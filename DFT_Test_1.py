# Import necessary PySCF modules
from pyscf import gto, dft, scf

# Define the molecule
mol = gto.M(
    atom='Ca 0.0 0.0 0.0; C 0.0 0.0 3.04; O 0.0 2.39 4.43; O 0.0 -2.39 4.43; O 0.0 0.0 -5.07',
    basis='def2-TZVPP',
    symmetry=True,
    verbose=4
)

# Define the DFT method
method = dft.RKS(mol)
method.xc = 'PBE'

# Compute the ground state energy
energy = method.kernel()

# Optimize the geometry of the molecule
scf.RHF(mol).kernel()

# Compute the vibrational frequencies
vib = mol.vibronic_spectrum(method, 10000)

# Print the results to the console
print('Ground state energy: {} a.u.'.format(energy))
print('Vibrational modes: {} cm^-1'.format(vib))
