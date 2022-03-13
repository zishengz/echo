# EChO

**E**lectro-**Ch**emical **O**ptimizer

`EChO` is an toolkit for local optimization and transition state search of electrochemical interfaces under a fixed electrode potential and in a grand canonical ensemble of electrons.


## Requirements

- `Python` 3.6 or higher
- `ASE` 3.22.1 or higher
- `VASP` 5.4.1 or higher, with `VASPsol` add-on


## Installation

Simply download and unzip the tarball of the code, and then add the directory to your `PYTHONPATH`:
```bash
wget https://github.com/zishengz/echo/archive/refs/heads/main.zip
unzip main.zip
export PYTHONPATH=$PYTHONPATH:`pwd`/echo-main
rm main.zip
```

In the jobs submission script, besides the `PYTHONPATH`, you should also set some [environment variables](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#environment-variables) for the `ASE` interface to `VASP`:
```bash
export PYTHONPATH=$PYTHONPATH:/xxx/echo-main
export VASPHOME=/xxx/vasp.x.x.x/bin
export VASP_COMMAND="mpirun -n $NSLOTS $VASPHOME/vasp_std"
export VASP_PP_PATH=/xxx/potentials/
```

The `VASP_COMMAND` would vary from machine to machine, depending on job scheduler and number of cores to use. You may need to adjust it yourself (and the script header, of course).


## Tutorials

Coming soon...
