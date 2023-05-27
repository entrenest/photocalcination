import numpy as np
from pyscf import gto, scf, grad, geomopt, tools, qmmm, cc, lib

# Define the Calcium Carbonate system
mol = gto.Mole()
mol.atom = '''
Ca  0.0000  0.0000  0.0000
C   1.6850  0.0000  0.0000
O   2.4150  0.9430  0.0000
O   2.4150 -0.9430  0.0000
C  -1.6850  0.0000  0.0000
O  -2.4150 -0.9430  0.0000
O  -2.4150  0.9430  0.0000
'''
mol.basis = '6-31g'
mol.verbose = 4
mol.build()

# Run a single point energy calculation
mf = scf.RHF(mol)
mf.kernel()

# Calculate the force constants for the system
grad_method = grad.RHF(mf)
fc = grad_method.hessian().reshape(-1, 3)

# Define the initial and final states
initial_state = np.array([0.0, 0.0, 0.0, 1.6850, 0.0, 0.0, 2.4150, 0.9430, 0.0, 2.4150, -0.9430, 0.0, -1.6850, 0.0, 0.0, -2.4150, -0.9430, 0.0, -2.4150, 0.9430, 0.0])
final_state = np.array([0.0, 0.0, 0.0, 2.4150, 0.9430, 0.0, 2.4150, -0.9430, 0.0, -2.4150, -0.9430, 0.0, -2.4150, 0.9430, 0.0, -1.6850, 0.0, 0.0, 1.6850, 0.0, 0.0])

# Define the intermediate states using interpolation
nimages = 5 # number of intermediate states
images = [initial_state] + [np.linspace(initial_state, final_state, nimages + 2)[i] for i in range(1, nimages + 1)] + [final_state]

# Define the NEB object and optimize the states
neb = tools.neb.NEB(fc, images)
neb.interpolate()
neb.optimizer = geomopt.BFGS(neb)
neb.kernel()

# Print out the optimized intermediate states
for i, image in enumerate(neb.images):
    print(f"Image {i+1} Energy: {neb.energies[i]:.6f}")
    print(image.reshape(-1, 3))

