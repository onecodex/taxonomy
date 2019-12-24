
build:
	maturin build --cargo-extra-args="--features=python" --release --strip
	docker run --rm -v $(pwd):/io konstin2/maturin:master build --cargo-extra-args="--features=python" --release --strip
test:
	cargo test --features python_test
	cargo clippy --features python_test
	cargo fmt -- --check
coverage:
	docker run --security-opt seccomp=unconfined -v "${PWD}:/volume" xd009642/tarpaulin bash -c "apt-get update -y && apt-get install python3-all-dev -y && cargo tarpaulin --features python_test -v"
