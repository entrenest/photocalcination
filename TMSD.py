import numpy as np
import pyscf

# Define the molecule
mol = pyscf.gto.Mole()
mol.atom = '''
    C   0.0000000  0.0000000  0.0000000
    O   1.1628800  0.0000000  0.0000000
    O  -1.1628800  0.0000000  0.0000000
    '''
mol.basis = '6-31g'
mol.build()

# Define the calculation parameters
mf = pyscf.scf.RHF(mol)
mf.conv_tol = 1e-9
mf.kernel()

# Set up the TMSD calculation
tmsd = pyscf.tddft.TDA(mf)

# Define the number of excited states to calculate
nstates = 10

# Calculate the TMSD for each excited state
tmsd_results = []
for state in range(nstates):
    # Calculate the excited state
    tmsd.exciton = state
    tmsd.kernel()
    # Calculate the TMSD for the excited state
    tmsd_res = tmsd.tmsd()
    tmsd_results.append(tmsd_res)

# Print the TMSD results
print('TMSD Results:')
for i, tmsd_res in enumerate(tmsd_results):
    print(f'State {i+1}: {tmsd_res:.6f} bohr^2')

