import numpy as np
from pyscf import gto, scf, dft, df, lib
from pyscf.prop.mulliken import MullikenPopulation
from pyscf.prop.nmr import NMR
from pyscf.tdscf import rhf

# Define the system: Calcium Carbonate
mol = gto.M(atom='''
Ca 0.0 0.0 0.0
C 1.781 1.781 1.781
O -1.781 -1.781 1.781
O -1.781 1.781 -1.781
O 1.781 -1.781 -1.781''', basis='cc-pvdz')

# Set up the molecular dynamics simulation
timestep = 0.5 # in femtoseconds
equilibration_time = 10 # in picoseconds
simulation_time = 50 # in picoseconds
temperature = 300 # in Kelvin

# Equilibrate the system
mf = dft.RKS(mol)
mf.xc = 'pbe,pbe'
mf.kernel()

# Set up the perturbation for thermal conductivity calculation
natom = mol.natm
nstep = int(simulation_time / timestep)
dt = timestep / lib.param.BOHR
vel = np.random.normal(size=(natom,3)) * np.sqrt(temperature / mf.mol.atom_mass_list.reshape(-1,1)) / lib.param.BOHR
vel -= np.average(vel, axis=0)
mol.set_velocities(vel)
mol.set_magnetic_field([0,0,0.01])
mol.set_charge([0,]*natom)
mol.set_spin([0,]*natom)
mol.set_initial_magnetic_moments()
mol.set_common_origin([0,0,0])

# Run the molecular dynamics simulation for thermal conductivity calculation
thermo = []
for istep in range(nstep):
    mol.update(velocity_verlet, dt)
    if istep * timestep >= equilibration_time:
        thermo.append(mol.thermodynamic_properties())

# Calculate the thermal conductivity using the Green-Kubo formula
kJ_to_eV = 0.010364 # Conversion factor from kilojoules per mole to electronvolts
heat_flux = np.array([prop.heatflux for prop in thermo])
corr_func = np.correlate(heat_flux[:,2], heat_flux[:,2], mode='full')
corr_time = np.arange(len(corr_func)) * timestep / 1000 # in picoseconds
integrand = corr_func[len(corr_func)//2:] / mol.vol * kJ_to_eV**2
lambda_thermal = np.trapz(integrand, corr_time[len(corr_time)//2:])

# Set up the perturbation for diffusion coefficient calculation
nstep = int(simulation_time / timestep)
dt = timestep / lib.param.BOHR
vel = np.random.normal(size=(natom,3)) * np.sqrt(temperature / mf.mol.atom_mass_list.reshape(-1,1)) / lib.param.BOHR
vel -= np.average(vel, axis=0)
mol.set_velocities(vel)
mol.set_magnetic_field([0,0,0])
mol.set_charge([0,]*natom)
mol.set_spin([0,]*natom)
mol.set_initial_magnetic_moments()
mol.set_common_origin([0,0,0])

# Run the molecular dynamics simulation for diffusion coefficient calculation
msd = np.zeros((nstep,))
for istep in range(n

