import os
from ase import Atoms
from mattersim.forcefield import MatterSimCalculator

def compute_isolated_energy(symbol):
    atom = Atoms(symbol, positions=[[0, 0, 0]], cell=[15, 15, 15], pbc=[False, False, False])
    atom.calc = MatterSimCalculator()
    return atom.get_potential_energy()

print("Calculating reference energies...")

energy_S = compute_isolated_energy("S")
energy_Se = compute_isolated_energy("Se")

print("Energy of isolated S (eV):", energy_S)
print("Energy of isolated Se (eV):", energy_Se)