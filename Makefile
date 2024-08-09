.PHONY: all
all: check test

.PHONY: test
bench:
	cargo +nightly bench --no-default-features --features no_std
	cargo +nightly bench --no-default-features --features std

.PHONY: check
check:
	cargo clippy --no-default-features --features no_std -- -D warnings
	cargo clippy --no-default-features --features std -- -D warnings
	cargo fmt --all -- --check

.PHONY: test
test:
	cargo test --no-default-features --features no_std
	cargo test --no-default-features --features std
