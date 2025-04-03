# EChO tutorial
Please firstly make sure you have installed properly the `EChO` and `VASPsol` and have set the related environment variables.


## Fixed-potential geometry optimization
The related files are in `examples/geomopt_fixpot`.

### Inputs
You will need to:
- Provide the initial geometry in VASP format (or whatever format that includes cell information and constraints). In the example they are `IS.vasp` and `FS.vasp`.
- Set VASP parameters in `geomopt_fixpot.py`.
- Set `pot_target` to the desired potential values, with `u_ref` being the SHE reference level.

### Outputs
Supposing we are running the job for `IS.vasp`:
- `IS_conv.db` contains the converged geometry with all calculated properties (potential energy, electronic free energy, forces, net charge, potential).
- `IS.log` stores the iteration printout of the geometry optimization.
- `IS.traj` stores the trajectory of the geometry optimization, with forces and potential energy.
- `LOG.xxx` is the output by job scheduler, which contains `EChO` printout. You may redirect it to any file in the submission script.

### Comments
- We recommend using initial geometry optimized under fixed-charge condition (regular approach by VASP) to speed up the force convergence, unless you would like to study a specific configuration that is only stabilized under applied potential.
- Recommended values of `u_ref` are 4.44 V (IUPAC value, theoretically self-consistent) and 4.6 V (benchmarked against exprimental PZC by Hennig groups).
- Multiple fixed-potential geometry optimization jobs can be run in the same directory since a working directory will be created for each job.
- The way `ASE` works with `VASP` is by writing/reading wavefunctions, so the job will be IO-intensive and takes up a lot of disk space. 

## Fixed-potential nudged elastic band calculation
The related files are in `examples/cineb_fixpot`.

### Inputs
You will need to:
- Provide `IS_conv.db` and `FS_conv.db` from **converged** fixed-potential geometry optimization of the initial and final state.
- You can provide a structure trajectory file (here we use `band_fixchg.db`) as the initial guess of the geometries along the band. Otherwise IDPP method will be used to interpolate between initial and final state geometries.
- Set `pot_target` to the desired potential values, with `u_ref` being the SHE reference level.

### Outputs
- `neb_conv.db` contains the converged band with all calculated properties (potential energy, electronic free energy, forces, net charge, potential).
- `band_init.db` stores the initial band.
- `neb.log` stores the iteration printout of the geometry optimization.
- `neb.traj` stores the trajectory of the geometry optimization, with forces and potential energy. 
- `restart.db` contains the band from the last CINEB iteration.
- `restart.nelect_neu` contains the net `NELECT` from the last CINEB iteration.
- `LOG.xxx` is the iteration printout by `EChO`.


### Comments
- It is highly recommended to start from a band from converged (CI)NEB run under fixed-charge condition.
- Parralelization over images is not supported.


## Fixed-potential molecular dynamics (MD)
The related files are in `examples/md_fixpot`

### Inputs
You would only need a structure file and an EChO script.

EChO will use the ASE calculator interface of VASP for the force calls. Hence, you need to define the VASP Calculator in the main py script.

You also need to provide the potentiostating parameters, see below:
- pot_target: the target potential, relative to the provided reference potential
- pot_step: a multiplier factor for potentiostating, in unit of e/V. It will be multiplied to the difference between current and target potentials to determine the amount of fractional electrons to add/remove. 
- pot_ref: the reference potential in vacuum scale. 4.6 V is the benchmarked value for VASPsol.
- potentiostat_nsteps: the interval of performing potentiostating. At $n$, EChO does potentiostating every $n$ MD step. 

All other parameters are the same as in a regular Langevin MD.

### Outputs

The raw output file (EChO.log) written by the EChO script will contain all the essential potentiostating information such as net charge, potential, and electronic grand canonical free energy.

Other MD-related information, such as temperature and potential energy, are written by ASE to a separate file (md-fp.log).

The full trajectory will be written to a ASE datafile on the fly (typically very large, therefore not provided in the example set).


## Fixed-potential constrained MD + thermodynamic integration
The related files are in `examples/slowgrowth_fixpot`

This can be used to obtain the fixed-potential free energy profile via the [Blue Moon Ensemble method](https://www.vasp.at/wiki/index.php/Blue_moon_ensemble).

In this study, we use the internal slow-growth functionality and SHAKE constraints of VASP for each MD step. `EChO` will handle the electronic degree of freedom (i.e., potentiostating) and keep track of all key quantities.


### Inputs
In this section, we are not using the ASE calculator interface of VASP, but directly calling VASP instead. Hence, the minimal input set includes a full VASP input set:
- `POSCAR`, `POTCAR`, `KPOINTS`
- `INCAR` which calls **one** step of slow-growth. Please do not add NELECT because it will be rewritten by EChO in the 1st line for every FPMD step.
- `ICONST` which describe the collective variable to be constrained and sampled.
- `echo_sg.py` which sets the FPMD-related parameters.

### Outputs
Please keep the raw output file written by the EChO script, as it contains all the essential information such as collective variable, temperature, potential, net charge, potential energy, electronic grand canonical free energy, and free energy gradient.

The full trajectory will be written to a ASE datafile on the fly (typically very large, therefore not provided in the example set).

For data extraction and analysis, you may refer to `examples/slowgrowth_fixpot/analysis.ipynb` for an exmaple.




