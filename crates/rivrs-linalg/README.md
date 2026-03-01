# rivrs-linalg

Numerical linear algebra for the [Rivrs](https://github.com/jcallaham/rivrs-linalg) framework.

Rivrs is a symbolic-numeric compiler for Rust that can build, transform, and execute computational graphs.
**Rivrs is currently a work in progress and the main compiler is not yet publicly available**.

To enable scientific computing applications, Rivrs as a framework includes Rust-native implementations of key numerical algorithms.
Since these can be used independently of the symbolic compiler, they are released as standalone crates.

`rivrs-linalg` contains numerical linear algebra functionality that builds on and extends [faer](https://crates.io/crates/faer).
It supports shared-memory parallelism using faer's `Par` type.

Currently `rivrs-linalg` is limited to a re-export of the [`rivrs-sparse`](https://crates.io/crates/rivrs-sparse) crate as a feature-gated module (`sparse`).

## Usage

```toml
[dependencies]
rivrs-linalg = "0.1"
```

The core functionality of the `rivrs-sparse` library is currently a sparse symmetric indefinite direct solver, a variant of Cholesky decomposition that factorizes a symmetric indefinite matrix $A$ into:

$$
A = L D L^T,
$$

where $D$ is block-diagonal and $L$ is lower-triangular.
Once this factorization is complete it can be reused for efficient linear system solves.

```rust,ignore
use rivrs_linalg::sparse::symmetric::SparseLDLT;
```

## Features

| Feature  | Default | Description |
|----------|---------|-------------|
| `sparse` | Yes     | Sparse linear algebra |

## License

Apache-2.0. See [LICENSE](LICENSE) for details.

For full documentation, see the [repository README](https://github.com/jcallaham/rivrs-linalg).
For solver documentation, see [`rivrs-sparse`](https://docs.rs/rivrs-sparse).
