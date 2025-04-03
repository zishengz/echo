from ase.calculators.vasp import Vasp
from ase.io import read
from echo.md import langevin_fixpot
import sys

vaspsol = Vasp(
    lcharg=False,
    npar=4,
    xc='pbe',
    ivdw=11,
    ismear=0,
    sigma=0.1,
    encut=300,
    lreal='Auto',
    prec='Low',
    # isym=0,
    algo='Fast',
    nelm=300,
    ediff=1e-5,
    gamma=True,
    kpts=[1, 1, 1],
    charge=0,
    lsol=True, eb_k=78.4, tau=0, lambda_d_k=3,
    setups='recommended',
    txt='vasp.out',
    directory='tmp')

geom_inp = read(sys.argv[1])

langevin_fixpot(
    atoms=geom_inp,
    timestep_fs=1,
    temperature_K=300,
    nsteps=1000,
    friction=0.1,
    label='md-fp',
    pot_target=0,
    pot_step=1,
    pot_ref=4.6,
    potentiostat_nsteps=5,
    calculator=vaspsol
)

