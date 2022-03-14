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

## Fixed-potential nedged elastic band calculation
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
- `LOG.xxx` is the ieration printout by `EChO`.


### Comments
- It is highly recommended to start from a band from converged (CI)NEB run under fixed-charge condition.
- Parralelization over images is not supported.

