from ase.calculators.vasp import Vasp
from ase.io import read, write
from ase.db import connect
from ase.optimize import FIRE
import sys

from echo.neb import neb_fix_pot

vaspsol = Vasp(
    lcharg=False,
    npar=8,
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
    setups='recommended',
    txt='vasp.out')

# Read from converged db file from fixpot geomopt
fn_ini = sys.argv[1]
with connect(fn_ini) as db:
    for r in db.select():
        geom_ini = r.toatoms()
        nelect_net_ini = r.nelect_net
        pot_ini = r.pot
        gcfe_ini = r.gcfe
fn_fin = sys.argv[2]
with connect(fn_fin) as db:
    for r in db.select():
        geom_fin = r.toatoms()
        nelect_net_fin = r.nelect_net
        pot_fin = r.pot
        gcfe_fin = r.gcfe

band_inp = None
# try:
#     n_images = len(read('band_init.db', index=':'))
#     band_inp = read('neb.traj', index=f'{-n_images}:')
# except:
#     band_inp = None
band_inp = read('band_fixchg.db', index=':')

band_conv = neb_fix_pot(
    geom_ini,
    geom_fin,
    pot_ini=pot_ini,
    pot_fin=pot_fin,
    nelect_net_ini=nelect_net_ini,
    nelect_net_fin=nelect_net_fin,
    gcfe_ini=gcfe_ini,
    gcfe_fin=gcfe_fin,
    pot_target=0,
    pot_ref=4.6,
    climb=True,
    calculator=vaspsol,
    optimizer=FIRE,
    fmax=0.02,
    pot_conv=0.001,
    band_inp=band_inp
)

write('geom_conv.db', band_conv)
