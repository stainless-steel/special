.PHONY: all
all: check test

.PHONY: test
bench:
	cargo +nightly bench
	cargo +nightly bench --no-default-features

.PHONY: check
check:
	cargo clippy -- -D warnings
	cargo clippy --no-default-features -- -D warnings
	cargo fmt --all -- --check

.PHONY: test
test:
	cargo test
	cargo test --no-default-features
