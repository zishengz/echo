from echo.data import pp_dict
from ase.db import connect
from ase.io import read, write
from ase.geometry import find_mic
from ase.calculators.vasp import Vasp
import numpy as np
import os



def get_fmax(forces):
    return np.sqrt((forces**2).sum(axis=1).max())


def get_nelect_neu(poscar='POSCAR', potcar='POTCAR', dirName='.'):
    elem_list = [l.split()[-2]
                 for l in open(f'{dirName}/{potcar}').readlines() if 'TITEL' in l]
    elem_list = [e.split('_')[0] if '_' in e else e for e in elem_list]
    zval_list = [eval(l.split()[-4])
                 for l in open(f'{dirName}/{potcar}').readlines() if 'ZVAL' in l]
    atom_list = read(f'{dirName}/{poscar}').get_chemical_symbols()
    return sum([zval_list[elem_list.index(a)] for a in atom_list])


def get_nelect_neu_ase(atoms):
    # Return neutral charge state via ASE
    # There must be a calculator attached to the Atoms!
    calc_tmp = Vasp()
    calc_tmp.fromdict(atoms.calc.todict())
    calc_tmp.set(directory='.', charge=0)
    calc_tmp.write_input(atoms)
    #nelect_neu = eval([l for l in open('INCAR').readlines() if 'NELECT' in l][0].split()[-1])
    nelect_neu = get_nelect_neu()
    for f in ['POSCAR', 'POTCAR', 'INCAR', 'KPOINTS', 'ase-sort.dat']:
        os.remove(f)
    return nelect_neu

def get_nelect_dict_pp(atoms, setting='recommended', pp_dir='potpaw_PBE'):
    list_elems_unique = list(set(atoms.get_chemical_symbols()))
    vasp_pp_path = os.environ['VASP_PP_PATH']
    nelect_dict = {}
    for elem in list_elems_unique:
        try:
            elem_pp = pp_dict[setting][elem]
        except:
            elem_pp = elem
        nelect_dict[elem] = eval([l for l in open(f'{vasp_pp_path}/{pp_dir}/{elem_pp}/POTCAR').readlines() if 'ZVAL' in l][0].split()[5])
    return nelect_dict


def get_nelect_neu_pp(atoms, setting='recommended', pp_dir='potpaw_PBE'):
    list_elems = atoms.get_chemical_symbols()
    nelect_dict = get_nelect_dict_pp(atoms=atoms, setting=setting, pp_dir=pp_dir)
    nelect_neu = 0
    for elem in list_elems:
        nelect_neu += nelect_dict[elem]
    return nelect_neu

def get_nelect_incar(workdir):
    return eval([l for l in open(f'{workdir}/INCAR').readlines() if 'NELECT' in l][0].split()[-1])

def extract_VASPsol(dirName='.', outname='vasp.out'):
    tmpHome = os.getcwd()
    os.chdir(dirName)
    fermi_shift = eval([l for l in open(outname).readlines()
                        if 'FERMI_SHIFT' in l][-1].split()[2])
    e_fermi = eval([l for l in open('OUTCAR').readlines()
                    if 'E-fermi' in l][-1].split()[2])
    nelect_now = eval([l for l in open('OUTCAR').readlines()
                       if 'NELECT' in l][-1].split()[2])
    energy_dft = eval([l for l in open('OSZICAR').readlines()
                       if 'F=' in l][-1].split()[4])
    os.chdir(tmpHome)
    return nelect_now, energy_dft, e_fermi, fermi_shift


def get_GCDFT(u_ref=4.6, poscar='POSCAR', potcar='POTCAR', dirName='.', outname='vasp.out'):
    nelect_neu = get_nelect_neu(poscar, potcar, dirName)
    nelect_now, energy_dft, e_fermi, fermi_shift = extract_VASPsol(
        dirName, outname)
    nelect_net = nelect_now - nelect_neu
    work_function = -e_fermi - fermi_shift
    u_she = work_function - u_ref
    energy_gcdft = energy_dft + nelect_net*(fermi_shift+work_function)
    return nelect_net, u_she, energy_gcdft


def get_GCDFT_neb(workdirs, u_ref=4.6, poscar='POSCAR', potcar='POTCAR', outname='vasp.out'):
    nelect_neu = get_nelect_neu(f'img001/{poscar}', f'img001/{potcar}')
    neb_nelect_net = []
    neb_u_she = []
    neb_energy_gcdft = []
    for dir in workdirs:
        nelect_now, energy_dft, e_fermi, fermi_shift = extract_VASPsol(
            dir, outname)
        nelect_net = nelect_now - nelect_neu
        work_function = -e_fermi - fermi_shift
        u_she = work_function - u_ref
        energy_gcdft = energy_dft + nelect_net*(fermi_shift+work_function)
        # append results
        neb_nelect_net.append(nelect_net)
        neb_u_she.append(u_she)
        neb_energy_gcdft.append(energy_gcdft)
    return neb_nelect_net, neb_u_she, neb_energy_gcdft


def get_rmsd(s, s_ref, list_ind=None):
    if list_ind is not None:
        tmp = s.positions[list_ind] - s_ref.positions[list_ind]
    else:
        tmp = s.positions - s_ref.positions
    if s.cell is not None and s.pbc is not None:
        tmp, _ = find_mic(tmp, s.cell, s.pbc)
    tmp **= 2
    return np.sqrt(tmp.sum(axis=1).mean())

def get_rssd(s, s_ref, list_ind=None):
    if list_ind is not None:
        tmp = s.positions[list_ind] - s_ref.positions[list_ind]
    else:
        tmp = s.positions - s_ref.positions
    if s.cell is not None and s.pbc is not None:
        tmp, _ = find_mic(tmp, s.cell, s.pbc)
    tmp **= 2
    return np.sqrt(tmp.sum(axis=1).sum())

# NEB band interval conversion


# NEB band transformation to diff IS/FS

