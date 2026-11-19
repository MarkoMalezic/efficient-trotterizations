import os

import numpy as np
from scipy.linalg import expm

from qiskit.circuit import QuantumCircuit
from methods import Heisenberg, Noise_analysis
from qiskit.circuit.library import PauliEvolutionGate


# Routines
def run_single(state, circuit, options, U_exact, noise=(1e-6, 1e-5), verbose=True):
  "Run a given circuit for a single state"

  results = Noise_analysis(state, circuit, options, U_exact, noise)
  if verbose:
    print("Trotter fidelity: Ft = ", results["trotter_fidelity"])
    print("Noisy fidelity: Fn = ", results["noisy_fidelity"])
    print("Noisy-Trotter fidelity: Fnt = ", results["noisy_trotter_fidelity"])
    print("Trotter error: error^2 = ", results["trotter_error2"])
    print("Depth: ", results["depth"])
    print("Size: ", results["size"])
    print("Operator counts: ", results["gates"])
  return results

def run_states(circuit, options, U_exact, noise=(1e-6, 1e-5), verbose=True):
  "Run a given circuit for all states in the computational basis"
  n = circuit.num_qubits
  states = [format(i, f"0{n}b") for i in range(2**n)]
  fidelities = np.zeros((len(states), 3))
  error2s = np.zeros((len(states), 1))
  for i, state in enumerate(states):
    print(f"Evolving state |{state}>, {i+1}/{2**n}", end="\r")
    results = Noise_analysis(state, circuit, options, U_exact, noise)
    fidelities[i, 0] = results["trotter_fidelity"]
    fidelities[i, 1] = results["noisy_fidelity"]
    fidelities[i, 2] = results["noisy_trotter_fidelity"]
    error2s[i, 0] = results["trotter_error2"]

  if verbose:
    print("Average Trotter fidelity: Ft = ", np.mean(fidelities[:, 0]))
    print("Average Noisy fidelity: Fn = ", np.mean(fidelities[:, 1]))
    print("Average Noisy-Trotter fidelity: Fnt = ", np.mean(fidelities[:, 2]))
    print("Average Trotter error: Et = ", np.sqrt(np.mean(error2s)))
    print("Depth: ", results["depth"])
    print("Size: ", results["size"])
    print("Operator counts: ", results["gates"])
  return fidelities, error2s, results["depth"], results["size"], results["gates"]


def save_step(circuit, options, reps, U_exact, save_path, noise=(1e-6, 1e-5)):
  "Save the results of running the circuit for a given number of repetitions"
  options["reps"] = reps
  fidelities, error2s, depth, size, gates = run_states(circuit, options, U_exact, noise=noise)
  # Save data
  os.makedirs(save_path, exist_ok=True)
  np.save(save_path + f"fidelities_Nt-{reps}.npy", fidelities)
  np.save(save_path + f"error2s_Nt-{reps}.npy", error2s)
  with open(save_path + f"gates_Nt-{reps}.txt", "w") as f:
    f.write(f"Depth: {depth}\n")
    f.write(f"Size: {size}\n")
    f.write(f"Gates: {gates}\n")


if __name__ == "__main__":

    # ---------------------- Model parameters ----------------------

    nq = 6   # number of spins/qubits

    # coupling constants
    Jx, Jy, Jz = 1.0, 0.0, 1.0
    # local field strength
    hz = np.array([0.05930859685756919, -0.06331304242132631, 0.05593819952253224, 0.019370032316011306, -0.010833448490313318, -0.08000501589123185])
    #hz = np.array([0.05930859685756919, -0.06331304242132631, 0.05593819952253224, 0.019370032316011306, -0.010833448490313318])
    # Boundary condition
    bound = "open" # "open" or "periodic"

    # 1) Build Hamiltonian
    H = Heisenberg(nq, Jx=Jx, Jy=Jy, Jz=Jz, hz=hz, bound=bound)

    # ----------------- Time evolution parameters ------------------

    t = 10.0              # total evolution time

    # Build the circuit for the time evolution operator U(t) = exp(-i t H)
    evo = PauliEvolutionGate(H, time=t)
    circuit = QuantumCircuit(nq)
    circuit.append(evo, list(np.arange(0, nq, 1)))

    # ---------------------- State evolution -----------------------

    # Exact time evolution using matrix exponentiation
    H_mat = H.to_matrix()
    U_exact = expm(-1j * t * H_mat)

    # ------------------------- Schemes ----------------------------

    # Custom schemes, given by order, cycles and ramp parameters c_vec
    n4_q6_recommended = {"cycles": 6, "order": 4, "c_vec": [0.07408257216024182, 0.2329230884208487, 0.29682056053401784, 0.1220869894365332, -0.35015363229365665, 0.1242404217420151], "reps": 100, "merge_single": True, "merge_steps": True}
    n6_q14_recommended = {"cycles": 14, "order": 6, "c_vec": [0.03725088315107759, 0.12060107932461332, 0.2660637966301682, 0.16366901594014166, 0.07131647756325317, 0.05811719881980566, 0.18870851618705248, -0.20001676266456045, 0.07414619113067222, 0.0873458895740801, 0.04423482019767251, -0.23082250310203417, -0.23719860775645235, 0.05658400500451008], "reps": 100, "merge_single": True, "merge_steps": True}
    n8_q21_recommended = {"cycles": 21, "order": 8, "c_vec": [0.041644811759980006, 0.08240131669062511, -0.10005367987707557, 0.056276640904364274, 0.05478929342750164, 0.052417975744997715, 0.20214313701451456, -0.15638705001185962, 0.04747687833330991, 0.07795810312185587, -0.21733485421642806, 0.07795810312185589, 0.04747687833330991, -0.15638705001185965, 0.20214313701451447, 0.052417975744997625, 0.05478929342750155, 0.0562766409043642, -0.10005367987707564, 0.08240131669062505, 0.041644811759979916], "reps": 100, "merge_single": True, "merge_steps": True}

    # Built-in schemes of qiskit-omelyan, selected by name
    leapfrog2 = {"name": "leapfrog2", "reps": 100, "merge_single": True, "merge_steps": True}
    omelyan2 = {"name": "omelyan2", "reps": 100, "merge_single": True, "merge_steps": True}
    forest_ruth4 = {"name": "forest_ruth4", "reps": 100, "merge_single": True, "merge_steps": True}
    malezic_ostmeyer4 = {"name": "malezic_ostmeyer4", "reps": 100, "merge_single": True, "merge_steps": True}
    malezic_ostmeyer4_minimal = {"name": "malezic_ostmeyer4_minimal", "reps": 100, "merge_single": True, "merge_steps": True}
    yoshida6 = {"name": "yoshida6", "reps": 100, "merge_single": True, "merge_steps": True}
    blanes_moan6 = {"name": "blanes_moan6", "reps": 100, "merge_single": True, "merge_steps": True}
    malezic_ostmeyer6 = {"name": "malezic_ostmeyer6", "reps": 100, "merge_single": True, "merge_steps": True}
    malezic_ostmeyer6_minimal = {"name": "malezic_ostmeyer6_minimal", "reps": 100, "merge_single": True, "merge_steps": True}
    morales8 = {"name": "morales8", "reps": 100, "merge_single": True, "merge_steps": True}
    malezic_ostmeyer8_q21 = {"name": "malezic_ostmeyer8_q21", "reps": 100, "merge_single": True, "merge_steps": True}
    malezic_ostmeyer8_q23 = {"name": "malezic_ostmeyer8_q23", "reps": 100, "merge_single": True, "merge_steps": True}
    malezic_ostmeyer8_q25 = {"name": "malezic_ostmeyer8_q25", "reps": 100, "merge_single": True, "merge_steps": True}
    malezic_ostmeyer8_q30 = {"name": "malezic_ostmeyer8_q30", "reps": 100, "merge_single": True, "merge_steps": True}

    # ------------------------- Noise model ------------------------

    eps = 1e-4
    noise = (0.1*eps, eps)   # depolarizing probabilities of (single-qubit, two-qubit) gates

    # ------------------------- Run --------------------------------

    # Scheme to run and its number of cycles q, which selects the list of Trotter steps in Nts/time_{t}/Nts_q{q}.dat
    scheme, q = malezic_ostmeyer8_q21, 21
    scheme_name = scheme.get("name", f"custom_n{scheme.get('order')}_q{scheme.get('cycles')}")

    # Paths relative to this script, so it can be run from any directory
    base_dir = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(base_dir, "Nts", f"time_{t}", f"Nts_q{q}.dat"), "r") as f:
      reps_list = [int(line) for line in f if line.strip()]
    # Separate folder per scheme, so runs of different schemes don't overwrite each other
    save_path = os.path.join(base_dir, "Data", f"time_{t}", scheme_name, "")

    #run_states(circuit, scheme, U_exact, noise=noise, verbose=True)
    for reps in reps_list:
      save_step(circuit, scheme, reps, U_exact, save_path=save_path, noise=noise)
