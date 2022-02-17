use std::io::Cursor;

use criterion::{criterion_group, criterion_main, Criterion};

use taxonomy::json;

fn prune_bench(c: &mut Criterion) {
    let build_json = include_str!("../tests/data/ncbi_subset_tax.json");
    let taxonomy = json::load(Cursor::new(build_json), None).expect("Error loading json");

    c.bench_function("prune", move |b| {
        b.iter(|| taxonomy.prune_to(&["65574"], false).expect("Error pruning"));
    });
}

criterion_group!(benches, prune_bench);
criterion_main!(benches);
