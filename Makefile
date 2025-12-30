features := elliptic lambert_w no_std std

.PHONY: all
all: check test

.PHONY: bench
bench:
	cargo +nightly bench --all-features

.PHONY: check
check: $(addprefix check-,$(features))
	cargo clippy --all-features -- -D warnings
	cargo fmt --all -- --check

.PHONY: $(addprefix check-,$(features))
$(addprefix check-,$(features)): check-%:
	cargo clippy --features $* -- -D warnings

.PHONY: test
test: $(addprefix test-,$(features))
	cargo test --all-features

.PHONY: $(addprefix test-,$(features))
$(addprefix test-,$(features)): test-%:
	cargo build --features $*
