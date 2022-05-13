use std::io::Cursor;

use criterion::{criterion_group, criterion_main, Criterion};

use taxonomy::{json::load, Taxonomy};

fn str_taxonomy(c: &mut Criterion) {
    let build_json = include_str!("../tests/data/ncbi_subset_tax.json");
    let taxonomy = load(Cursor::new(build_json), None).expect("Error loading json");

    c.bench_function("lineage str", move |b| {
        b.iter(|| {
            taxonomy
                .lca(
                    65574u32.to_string().as_str(),
                    160352u32.to_string().as_str(),
                )
                .unwrap()
                .parse::<u32>()
                .unwrap()
                .to_string()
        });
    });
}

fn u32_taxonomy(c: &mut Criterion) {
    let build_json = include_str!("../tests/data/ncbi_subset_tax.json");
    let taxonomy = load(Cursor::new(build_json), None).expect("Error loading json");

    c.bench_function("lineage u32", move |b| {
        b.iter(|| taxonomy.lca(1577, 828));
    });
}

criterion_group!(benches, str_taxonomy, u32_taxonomy);
criterion_main!(benches);
