# rivrs

Rivrs is a symbolic-numeric compiler for Rust that can build, transform, and execute computational graphs.

**Rivrs is currently a work in progress and the main compiler is not yet publicly available**.

Since the numerical implementations can be used independently of the compiler for scientific computing applications, these are released as standalone crates that are re-exported by `rivrs`.
Currently this is limited to [`rivrs-linalg`](https://crates.io/crates/rivrs-linalg), containing numerical linear algebra functionality that builds on and extends [faer](https://crates.io/crates/faer).

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
