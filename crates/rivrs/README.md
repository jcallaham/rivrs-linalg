# rivrs

Scientific computing framework for Rust.

This is an umbrella crate that re-exports domain libraries as feature-gated modules:

- **`linalg`** (default) — [`rivrs-linalg`](https://crates.io/crates/rivrs-linalg): Numerical linear algebra

Future modules: `ode`, `optimize`, `signal`, `mesh`.

## Usage

```toml
[dependencies]
rivrs = "0.1"
```

```rust,ignore
use rivrs::linalg::sparse::symmetric::SparseLDLT;
```

## Features

| Feature  | Default | Description |
|----------|---------|-------------|
| `linalg` | Yes     | Numerical linear algebra |

## License

Apache-2.0. See [LICENSE](LICENSE) for details.

For full documentation, see the [repository README](https://github.com/jcallaham/rivrs-linalg).
