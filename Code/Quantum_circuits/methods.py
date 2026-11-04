import numpy as np

from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, SparsePauliOp, state_fidelity
from qiskit.transpiler.passes.synthesis.high_level_synthesis import HLSConfig
from qiskit.transpiler import generate_preset_pass_manager
from qiskit_aer.noise import NoiseModel, depolarizing_error
from qiskit_aer import AerSimulator
from qiskit.quantum_info import DensityMatrix

# Model: Heisenberg XXZ chain with local field
def Heisenberg(n, Jx=1.0, Jy=1.0, Jz=1.0, hz=0.1, bound="open"):
  """Return Heisenberg-chain Hamiltonian as SparsePauliOp on n qubits (open or periodic boundary)."""

  hz = np.asarray(hz)

  if hz.ndim == 0:
    hz = np.full(n, hz)
  elif hz.shape != (n,):
    raise ValueError(f"hz must be a scalar or length-{n} array")
  
  paulis = []
  coeffs = []

  def two_body_term(p, i, j, c):
    # p is "XX", "YY", or "ZZ"
    s = ["I"] * n
    s[i] = p[0]
    s[j] = p[1]
    paulis.append("".join(s))
    coeffs.append(c)

  def one_body_term(p, i, c):
    # p is "X" or "Z"
    s = ["I"] * n
    s[i] = p
    paulis.append("".join(s))
    coeffs.append(c)

  # nearest-neighbor couplings
  for i in range(n - 1):
    if Jx != 0: two_body_term("XX", i, i + 1, Jx)
    if Jy != 0: two_body_term("YY", i, i + 1, Jy)
    if Jz != 0: two_body_term("ZZ", i, i + 1, Jz)
  
  if bound == "periodic" and n > 2:
    # periodic boundary condition
    if Jx != 0: two_body_term("XX", n - 1, 0, Jx)
    if Jy != 0: two_body_term("YY", n - 1, 0, Jy)
    if Jz != 0: two_body_term("ZZ", n - 1, 0, Jz)

  # local field
  for i in range(n):
    if hz[i] != 0: one_body_term("Z", i, hz[i])

  return SparsePauliOp(paulis, coeffs=np.asarray(coeffs, dtype=float))


# Time evolution error analysis
def Trotter_error(evo_circuit, U_exact):
  "Trotter error computed as sum of 2-norm errors over all basis states."
  # Generate all basis states
  n = evo_circuit.num_qubits
  states = [format(i, f"0{n}b") for i in range(2**n)]
  error2 = 0.0
  # Evolve each basis state with the Trotter circuit and compare to exact evolution
  for i, state in enumerate(states):
    print(f"Testing state |{state}>, {i+1}/{2**n}", end="\r")
    psi0 = Statevector.from_label(state)
    trotter_psi = psi0.evolve(evo_circuit)
    exact_psi = Statevector(U_exact @ psi0.data)
    error2 += np.linalg.norm(exact_psi.data - trotter_psi.data)**2
  # Normalize by number of states
  error = np.sqrt(error2 / 2**n)
  return error

def Sweep_steps(circuit, options, reps_list, U_exact):
  "Run the given circuit with different numbers of steps and compute the Trotter error for each."
  errors = np.array([])
  for reps in reps_list:
    print("\nNumber of steps: ", reps)
    options["reps"] = reps
    hls_config = HLSConfig(PauliEvolution=[("omelyan_trotter", options)])
    pass_manager = generate_preset_pass_manager(optimization_level=3, hls_config=hls_config)
    evo_circuit = pass_manager.run(circuit)
    error = Trotter_error(evo_circuit, U_exact)
    errors = np.append(errors, error)
  
  return errors


# Noise-only benchmarking support
def zero_angle_circuit(circuit):
  "Replace every gate by its zero-angle (identity) version, keeping the same gate sequence/names so the same per-gate noise still applies."
  zeroed = QuantumCircuit(*circuit.qregs)
  for instr in circuit.data:
    op = instr.operation
    if op.params:
      op = type(op)(*([0.0] * len(op.params)))
    zeroed.append(op, instr.qubits, instr.clbits)
  return zeroed


# Zero-noise extrapolation support
def fold_circuit(circuit, num_folds):
  "Locally fold every gate num_folds times (G -> G (G^-1 G)^num_folds), scaling the circuit's noise by (2*num_folds + 1)."
  folded = QuantumCircuit(*circuit.qregs)
  for instr in circuit.data:
    folded.append(instr.operation, instr.qubits, instr.clbits)
    for _ in range(num_folds):
      folded.append(instr.operation.inverse(), instr.qubits, instr.clbits)
      folded.append(instr.operation, instr.qubits, instr.clbits)
  return folded


# Time evolution noise analysis
def Noise_analysis(state, circuit, options, U_exact, noise=(1e-6, 1e-5), mitigate=False, num_folds=(0, 1, 2), noise_only=False):
  """Run the given circuit with and without noise, and compute the fidelity between the noisy state and the exact state.

  If mitigate=True, additionally zero-noise-extrapolate the noisy fidelity using global unitary folding
  (see fold_circuit) at the given num_folds (scale factors 2*num_folds+1) and mitiq's Richardson extrapolation.
  mitiq's own circuit-folding/conversion machinery is not used here: it round-trips Qiskit circuits through
  QASM2/cirq, which doesn't support the rxx/ryy/rzz gates this circuit is built from and falls back to a
  qiskit.transpile() call that crashes under mitiq's pinned numpy<2.0 (qiskit's Operator construction relies
  on numpy>=2.0 `copy=None` semantics). Folding the native gates ourselves sidesteps that entirely.

  If noise_only=True, additionally run a version of the circuit with every rotation angle zeroed out (see
  zero_angle_circuit), keeping the same gate sequence/names so the noise model still applies per gate. The
  ideal (noiseless) evolution of that circuit is the identity, so its fidelity to psi0 isolates the fidelity
  loss from accumulated gate noise alone, decoupled from the Trotter truncation error.
  """

  # Compute the exact state after evolution
  psi0 = Statevector.from_label(state)
  exact_psi = Statevector(U_exact @ psi0.data)

  # Build Trotterized circuit with given options
  hls_config = HLSConfig(PauliEvolution=[("omelyan_trotter", options)])
  pass_manager = generate_preset_pass_manager(optimization_level=0, hls_config=hls_config)
  evo_circuit = pass_manager.run(circuit)
  # Evolve the initial state without noise
  trotter_psi = psi0.evolve(evo_circuit)
  # Compute the fidelity between the Trotterized state and the exact state
  trotter_fid = state_fidelity(trotter_psi, exact_psi)
  trotter_error2 = np.linalg.norm(exact_psi.data - trotter_psi.data)**2
  #print(f"Fidelity between trotter and exact state: {trotter_fid:.6e}")

  # Build the noise model
  noise_model = NoiseModel()
  noise_model.add_all_qubit_quantum_error(depolarizing_error(noise[0], 1), ["rz"])
  noise_model.add_all_qubit_quantum_error(depolarizing_error(noise[1], 2), ["rxx", "ryy", "rzz"])
  # Simulate the circuit with noise
  sim = AerSimulator(method="density_matrix", noise_model=noise_model)

  def noisy_fidelity(evolved_circuit, reference):
    "Prepend the initial state and run evolved_circuit through the noisy simulator, returning the fidelity to reference."
    noisy_circuit = QuantumCircuit(evolved_circuit.num_qubits)
    for i, bit in enumerate(reversed(state)):
      if bit == "1":
        noisy_circuit.x(i)
    noisy_circuit.compose(evolved_circuit, inplace=True)
    noisy_circuit.save_density_matrix()
    result = sim.run(noisy_circuit).result()
    rho_noisy = DensityMatrix(result.data(0)["density_matrix"])
    return rho_noisy, state_fidelity(rho_noisy, reference)

  if noise_only:
    # The zero-angle circuit's ideal (noiseless) output is psi0 itself, not exact_psi:
    # compare against psi0 so the fidelity isolates noise, not the real evolution's distance from psi0.
    evo_circuit = zero_angle_circuit(evo_circuit)
    reference = psi0
  else:
    reference = exact_psi

  rho_noisy, noisy_fid = noisy_fidelity(evo_circuit, reference)
  #print(f"Fidelity between noisy and exact state: {noisy_fid:.6e}")
  noisy_trotter_fid = state_fidelity(rho_noisy, trotter_psi if not noise_only else psi0)
  #print(f"Fidelity between noisy and trotter state: {noisy_trotter_fid:.6e}")

  results = {
    "trotter_fidelity": trotter_fid,
    "noisy_fidelity": noisy_fid,
    "noisy_trotter_fidelity": noisy_trotter_fid,
    "trotter_error2": trotter_error2,
    "depth": evo_circuit.depth(),
    "size": evo_circuit.size(),
    "gates": evo_circuit.count_ops(),
  }

  if mitigate:
    from mitiq.zne.inference import RichardsonFactory  # only needed for mitigate=True; requires the .venv_mitiq environment
    factory = RichardsonFactory(scale_factors=[2 * k + 1 for k in num_folds])
    for k in num_folds:
      fid = noisy_fid if k == 0 else noisy_fidelity(fold_circuit(evo_circuit, k), reference)[1]
      factory.push({"scale_factor": 2 * k + 1}, fid)
    results["mitigated_fidelity"] = factory.reduce()

  return results

