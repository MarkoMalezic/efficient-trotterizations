# Noisy quantum circuits

This folder contains the scripts used to compare Trotter-Suzuki schemes as quantum circuits under gate noise, as presented in "Optimised 8th-Order Trotter-Suzuki Schemes for Quantum Simulation" ([arXiv:2611.27896](https://arxiv.org/abs/2611.27896)).
The circuits are synthesised with our Qiskit plugin [qiskit-omelyan](https://github.com/MarkoMalezic/qiskit-omelyan) ([10.5281/zenodo.24371684](https://doi.org/10.5281/zenodo.24371684)), based on the algorithm described in our Lattice 2025 proceedings ([arXiv:2602.21145](https://arxiv.org/abs/2602.21145)). It implements arbitrary symmetric Trotter-Suzuki schemes (in the ramp notation with parameters c, see `README-Schemes.md`) as well as a collection of schemes from the literature and from our papers.
The circuits are simulated with Qiskit Aer, both without and with noise.

## Requirements

- Python >= 3.10
- `numpy`, `scipy`
- `qiskit` (>= 2.0), `qiskit-aer`
- `qiskit-omelyan` >= 1.0.0, installed from PyPI (earlier versions do not include the order 8 schemes and the `_minimal` schemes used here):
  ```bash
  pip install numpy scipy qiskit qiskit-aer
  pip install "qiskit-omelyan>=1.0.0"
  ```
- Optional: `mitiq`, only for zero-noise extrapolation (`mitigate=True` in `Noise_analysis`). Note that `mitiq` pins `numpy<2`.

## Usage

```bash
python noisy_circuits.py
```

The script has no command-line arguments. The model, the scheme and the noise strength are set in the `__main__` block of `noisy_circuits.py`:

- **Model**: the Heisenberg chain on `nq = 6` qubits
  $$H = \sum_{i} \left(J_x X_i X_{i+1} + J_y Y_i Y_{i+1} + J_z Z_i Z_{i+1}\right) + \sum_i h_i Z_i,$$
  with $J_x = J_z = 1$, $J_y = 0$, a fixed set of random local fields $h_i$ and open boundary conditions, evolved to time `t = 10`.
  The exact propagator $U(t) = e^{-iHt}$ is obtained by exponentiating the Hamiltonian matrix.
- **Scheme**: `scheme, q = malezic_ostmeyer8_q21, 21` selects the options passed to the plugin and the number of cycles q.
  The options are dictionaries with either the `name` of a built-in scheme of qiskit-omelyan, or `order`, `cycles` and `c_vec` for a custom scheme (see `n4_q6_recommended` etc.), together with the number of Trotter steps `reps` and the circuit optimisations `merge_single` and `merge_steps`.
- **Noise**: a depolarising channel after every gate, with probability `noise[0]` for single-qubit gates (`rz`) and `noise[1]` for two-qubit gates (`rxx`, `ryy`, `rzz`). The default is `noise = (1e-5, 1e-4)`.

For every number of Trotter steps $N_t$ listed in `Nts/time_{t}/Nts_q{q}.dat`, the script evolves all $2^{n_q}$ computational basis states and saves the results to `Data/time_{t}/<scheme name>/`:

| File | Content |
|---|---|
| `fidelities_Nt-{Nt}.npy` | Array of shape $(2^{n_q}, 3)$ with the fidelities per initial state: Trotter vs. exact ($F_t$), noisy vs. exact ($F_n$), noisy vs. Trotter ($F_{nt}$) |
| `error2s_Nt-{Nt}.npy` | Array of shape $(2^{n_q}, 1)$ with the squared Trotter error $\lVert \psi_\text{exact} - \psi_\text{Trotter} \rVert^2$ per initial state |
| `gates_Nt-{Nt}.txt` | Depth, size and gate counts of the synthesised circuit |

The noisy simulation uses density matrices, so the cost grows quickly with the number of qubits.
The simulation data is not included in the repository, but it will be gladly provided upon request.

## Files

| File | Content |
|---|---|
| `noisy_circuits.py` | Main script: model and scheme setup, and the routines `run_single` (one initial state), `run_states` (all basis states) and `save_step` (run all states for a given $N_t$ and save the results) |
| `methods.py` | `Heisenberg` (Hamiltonian as a `SparsePauliOp`), `Noise_analysis` (noiseless and noisy simulation of one initial state), and the helpers `Trotter_error` / `Sweep_steps` (noiseless Trotter error over all basis states for a list of $N_t$) |
| `Nts/time_10.0/Nts_q{q}.dat` | The numbers of Trotter steps $N_t$ (one per line) used for the scheme with q cycles at t = 10 |

`Noise_analysis` has two optional modes:
- `mitigate=True`: additionally extrapolates the noisy fidelity to zero noise (Richardson extrapolation from `mitiq`), with the noise scaled by folding every gate $G \to G (G^{-1} G)^k$.
- `noise_only=True`: runs the same circuit with all rotation angles set to zero, so that the fidelity loss from the gate noise is isolated from the Trotter error.
