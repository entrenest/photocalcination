# Import necessary PySCF and pymatgen modules
from pyscf import gto, dft, scf
from pyscf.scf import addons
from pymatgen.io.cif import CifParser

# Parse the CIF data
parser = CifParser("/home/phill/Desktop/Scripts/Quantum/CaCO3.cif")
structure = parser.get_structures()[0]

# Define the molecule
mol = gto.Mole()
mol.atom = [[str(site.specie), site.coords.tolist()] for site in structure]
mol.basis = 'cc-pvtz'
mol.symmetry = True
mol.verbose = 4
mol.build()

# Define the DFT method
method = dft.RKS(mol).newton()  # Use the Newton-Raphson method
method.xc = 'PBE'
method.max_cycle = 1000

# Compute the ground state energy
energy = method.kernel()

# Optimize the geometry of the molecule
opt = scf.RHF(mol).newton().run()

# Compute the vibrational frequencies
from pyscf.hessian import rhf as hrhf
hessian = hrhf.Hessian(opt)
freqs = hessian.kernel()

# Print the results to the console
print('Ground state energy: {} a.u.'.format(energy))
print('Vibrational frequencies: {} cm^-1'.format(freqs))
