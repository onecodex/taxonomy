# Taxonomy: Functional Tree Operations — API Design

## Context

Users need to store arbitrary data alongside taxonomy nodes and perform aggregation and transformation operations across the tree. This is motivated by use cases like computing subtree read counts in metagenomic analysis.

All operations are implemented in Rust and exposed to Python via PyO3. Operations on the tree return new `Taxonomy` objects (immutable/functional style), consistent with the existing `prune` method. `TaxonomyNode` objects remain independent value objects with no back-reference to the tree.

______________________________________________________________________

## Summary

| Operation | Traversal | Lambda | Complexity |
|---|---|---|---|
| `reduce_up` | post-order (leaves → root) | `f(node, [child_results]) -> result` | O(n) |
| `map_down` | pre-order (root → leaves) | `f(parent_result, node) -> result` | O(n) |

n = number of nodes in the subtree rooted at `node_id`.

______________________________________________________________________

## Data Access

### Reading — existing API

`TaxonomyNode` already exposes extra data via `__getitem__`. Data is populated from the underlying `data: Vec<HashMap<String, Value>>` field when a node is constructed.

```python
node = tax["562"]
node["readcount"]        # raises KeyError if key absent
node.get("readcount", 0) # returns default if absent — NEW
node.data                # full data dict — NEW
```

**New methods needed on `TaxonomyNode`:**

| Method | Complexity | Notes |
|---|---|---|
| `node.get(key, default=None)` | O(1) | safe read with fallback |
| `node.data` | O(d) | returns copy of data as Python dict, d = number of keys |

**Note:** `TaxonomyNode` is a snapshot — it reflects the tree state at the time it was constructed. Calling `set_data` after fetching a node does not update existing node references.

______________________________________________________________________

### Writing — new API

```python
tax.set_data(node_id: str, key: str, value) -> None
```

- Mutates the taxonomy in-place (consistent with `add_node`, `edit_node`)
- **O(1)**: hash map lookup by `node_id`, hash map insert for `key`

```python
tax.set_data("562", "readcount", 5)
tax["562"]["readcount"]  # 5
```

______________________________________________________________________

## Aggregation

### `reduce_up` — Aggregate from leaves to root

```python
tax.reduce_up(node_id: str, output_key: str, fn: Callable[[TaxonomyNode, List], result]) -> Taxonomy
```

- **O(n)** — visits every node in the subtree exactly once
- **Post-order** traversal: leaves visited before parents
- `fn(node, child_results) -> result`
  - `node`: the current `TaxonomyNode`
  - `child_results`: list of already-computed results from direct children (empty list for leaves)
- Stores result at **every node** under `output_key`
- Returns a **new Taxonomy** (original unchanged)
- No `initial` value — leaves handle the base case via `child_results == []`
- Chainable: results stored by one `reduce_up` are visible on nodes in the next

Mirrors `functools.reduce` conceptually: reduces the tree bottom-up.

```python
# Compute inclusive (clade) read counts — equivalent to Kraken's "clade_reads"
annotated = tax.reduce_up("1", "clade_reads",
    lambda node, child_results: node.get("readcount", 0) + sum(child_results))
annotated["562"]["clade_reads"]   # all reads in the E. coli clade
annotated["1224"]["clade_reads"]  # all reads in Proteobacteria

# Count detected species per clade
tax.reduce_up("1", "detected_species",
    lambda node, child_results: sum(child_results) + (1 if node.rank == "species" and node.get("readcount", 0) > 0 else 0))

# Compute relative abundance (chained)
annotated = tax.reduce_up("1", "clade_reads",
    lambda node, child_results: node.get("readcount", 0) + sum(child_results))
annotated.reduce_up("1", "relative_abundance",
    lambda node, child_results: node.get("readcount", 0) / annotated["1"]["clade_reads"])
```

______________________________________________________________________

### `map_down` — Propagate values from root to leaves

```python
tax.map_down(node_id: str, output_key: str, initial, fn: Callable[[parent_result, TaxonomyNode], result]) -> Taxonomy
```

- **O(n)** — visits every node in the subtree exactly once
- **Pre-order** traversal: parents visited before children
- `fn(parent_result, node) -> result`
  - `parent_result`: result stored at the parent (or `initial` for the root node)
  - `node`: the current `TaxonomyNode`
- Stores result at **every node** under `output_key`
- Returns a **new Taxonomy**
- Chainable with `reduce_up` and `map_down`

Mirrors Python's `map` conceptually: transforms each node using context flowing from its parent.

```python
# Build full lineage string for every node (QIIME-style taxonomy strings)
tax.map_down("1", "lineage", "",
    lambda parent_lineage, node: f"{parent_lineage};{node.name}" if parent_lineage else node.name)
# tax["562"]["lineage"]
# → "Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacterales;Enterobacteriaceae;Escherichia;Escherichia coli"

# Compute depth of every node
tax.map_down("1", "depth", 0,
    lambda parent_depth, node: parent_depth + 1)

# Propagate cumulative branch length from root
tax.map_down("1", "distance_from_root", 0.0,
    lambda parent_dist, node: parent_dist + node["branch_length"])
```

______________________________________________________________________

## Performance Notes

The lambda receives a full `TaxonomyNode` on every call, which currently requires allocating and populating a Python object per node (string copies for `id`, `name`, `rank`, `parent`, plus all data keys). For large trees (e.g. NCBI ~2M nodes) this has meaningful overhead. Two future optimization paths:

- **Zero-copy node**: pass a borrowed view backed by a pointer into the Rust tree (safe during traversal since the tree is not mutated), avoiding all allocations
- **Built-in Rust-native ops** (`sum`, `count`, `max`, `min`): bypass the lambda entirely for common cases

Both are deferred until the API is validated.

______________________________________________________________________

## Comparison to NetworkX and ete3

This library was written as a replacement for NetworkX for taxonomy use cases. Neither NetworkX nor ete3 have built-in equivalents of `reduce_up` or `map_down` — both require manual traversal loops.

### `reduce_up`

**NetworkX:**

```python
def reduce_up(G, root, fn):
    for node_id in nx.dfs_postorder_nodes(G, root):
        child_results = [G.nodes[c]["_result"] for c in G.successors(node_id)]
        G.nodes[node_id]["_result"] = fn(G.nodes[node_id], child_results)
```

**ete3:**

```python
for node in tree.traverse("postorder"):
    child_results = [c.clade_reads for c in node.children]
    node.clade_reads = node.readcount + sum(child_results)
```

**taxonomy:**

```python
annotated = tax.reduce_up("1", "clade_reads",
    lambda node, child_results: node.get("readcount", 0) + sum(child_results))
```

______________________________________________________________________

### `map_down`

**NetworkX:**

```python
def map_down(G, root, initial, fn):
    for node_id in nx.dfs_preorder_nodes(G, root):
        parents = list(G.predecessors(node_id))
        parent_result = G.nodes[parents[0]]["_result"] if parents else initial
        G.nodes[node_id]["_result"] = fn(parent_result, G.nodes[node_id])
```

**ete3:**

```python
for node in tree.traverse("preorder"):
    parent_lineage = node.up.lineage if not node.is_root() else ""
    node.lineage = f"{parent_lineage};{node.name}" if parent_lineage else node.name
```

**taxonomy:**

```python
annotated = tax.map_down("1", "lineage", "",
    lambda parent_lineage, node: f"{parent_lineage};{node.name}" if parent_lineage else node.name)
```

______________________________________________________________________

Key differences from NetworkX:

- NetworkX uses `DiGraph` with dict-style node attributes; this library uses typed `TaxonomyNode` objects with rank, name, and parent built in
- NetworkX has no concept of taxonomic rank, lineage, or LCA — these require manual implementation
- This library is implemented in Rust; NetworkX is pure Python

Key differences from ete3:

- ete3 uses attribute access (`node.readcount`); this library uses `node["readcount"]`
- ete3 is pure Python; this library is implemented in Rust
- ete3 has richer phylogenetic features (branch support, evolutionary models); this library is optimized for large taxonomic trees (NCBI ~2M nodes)

______________________________________________________________________

## Deferred

- Built-in Rust-native `sum`, `count`, `max`, `min` (optimization, post-validation)
- `map(output_key, fn)` — transform data values per node without aggregation
