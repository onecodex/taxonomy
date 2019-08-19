use std::io::Cursor;

use criterion::{criterion_group, criterion_main, Criterion};

use taxonomy::edit::prune_to;
use taxonomy::formats::json::load_json;

fn prune_bench(c: &mut Criterion) {
    let build_json = include_str!("../tests/data/ncbi_subset_tax.json");
    let taxonomy = load_json(Cursor::new(build_json), None).expect("Error loading json");
    let tax_id = taxonomy.to_internal_id("65574").expect("Tax ID not found");

    c.bench_function("prune", move |b| {
        b.iter(|| prune_to(&taxonomy, &[tax_id], false).expect("Error pruning"));
    });
}

criterion_group!(benches, prune_bench);
criterion_main!(benches);
