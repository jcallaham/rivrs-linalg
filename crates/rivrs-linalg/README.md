# rivrs-linalg

Numerical linear algebra for the [Rivrs](https://github.com/jcallaham/rivrs-linalg) scientific computing framework.

This is a facade crate that re-exports domain crates as feature-gated modules:

- **`sparse`** (default) — [`rivrs-sparse`](https://crates.io/crates/rivrs-sparse): Sparse symmetric indefinite solver (SSIDS)

## Usage

```toml
[dependencies]
rivrs-linalg = "0.1"
```

```rust,ignore
use rivrs_linalg::sparse::symmetric::SparseLDLT;
```

## Features

| Feature  | Default | Description |
|----------|---------|-------------|
| `sparse` | Yes     | Sparse symmetric indefinite solver |

## License

Apache-2.0. See [LICENSE](LICENSE) for details.

For full documentation, see the [repository README](https://github.com/jcallaham/rivrs-linalg).
For solver documentation, see [`rivrs-sparse`](https://docs.rs/rivrs-sparse).
