from echo.vaspmd import slow_growth
from ase.io import read, write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution


geom_inp = read('POSCAR')

slow_growth(
    geom_inp,
    nsteps = 6000,
    pot_target = -0.0591593 * 1,
    pot_ref = 4.6,
    pot_step = 928 * 0.001,
)


