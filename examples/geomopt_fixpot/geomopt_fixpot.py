
from ase.calculators.vasp import Vasp
from ase.io import read, write
from ase.optimize import LBFGSLineSearch
import sys

from echo.opt import opt_fix_pot

vaspsol = Vasp(
    lcharg=False,
    npar=4,
    xc='rpbe',
    # ivdw=11,
    ismear=0,
    sigma=0.1,
    encut=400,
    lreal='Auto',
    prec='Normal',
    # isym=0,
    algo='Fast',
    nelm=300,
    ediff=1e-5,
    gamma=True,
    kpts=[1, 1, 1],
    charge=0,
    lsol=True, eb_k=78.4, tau=0, lambda_d_k=3,
    setups='recommended')

fn_inp = sys.argv[1]
geom_inp = read(fn_inp)
label_inp = fn_inp.split('.')[0]

geom_out = opt_fix_pot(
    geom_inp,
    label=label_inp,
    pot_target=0.0,
    pot_ref=4.6,
#    pot_step=1.5,
    calculator=vaspsol,
    optimizer=LBFGSLineSearch,
    fmax=0.02,
    pot_conv=0.001
)

write(f'{label_inp}-opt.vasp', geom_out)
