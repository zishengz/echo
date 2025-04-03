# EChO

**E**lectro-**Ch**emical **O**ptimizer

`EChO` is an toolkit for local optimization and transition state search of electrochemical interfaces under a fixed electrode potential and in a grand canonical ensemble of electrons.

Please make sure to cite the following papers if you use `EChO` in your research:

- J. Am. Chem. Soc. 2024, 146, 14, 9623–9630. [[link]](https://doi.org/10.1021/jacs.3c12934)  [[bibtex]](https://scholar.googleusercontent.com/scholar.bib?q=info:qFJ9T1LlKo4J:scholar.google.com/&output=citation&scisdr=ClGobQuFEJXb1FnTVG0:AFWwaeYAAAAAZ-3VTG2hHFvnWSI4mdgEc7_TCo8&scisig=AFWwaeYAAAAAZ-3VTAgCu-zQ_y9aHxgAZxvbrV4&scisf=4&ct=citation&cd=-1&hl=en)

- J. Am. Chem. Soc. 2025, accepted. [[link]](10.1021/jacs.5c00775)  [[bibtex]](https://scholar.googleusercontent.com/scholar.bib?q=info:qFJ9T1LlKo4J:scholar.google.com/&output=citation&scisdr=ClGobQuFEJXb1FnTVG0:AFWwaeYAAAAAZ-3VTG2hHFvnWSI4mdgEc7_TCo8&scisig=AFWwaeYAAAAAZ-3VTAgCu-zQ_y9aHxgAZxvbrV4&scisf=4&ct=citation&cd=-1&hl=en)



## Requirements

- `Python` 3.6 or higher
- `ASE` 3.22.1 or higher
- `VASP` 5.4.1 or higher, with `VASPsol` add-on


## Installation

If your machine has Git installed, simply clone the repo to your local directory by:
```
git clone https://github.com/zishengz/echo
```
Or, you can also download and unzip the source code:
```
wget https://github.com/zishengz/echo/archive/refs/heads/main.zip
unzip main.zip
rm main.zip
mv main echo
```
After fetching the gocia repo, add it to your PYTHONPATH by:
```
export PYTHONPATH=$PYTHONPATH:`pwd`/echo/
```
Remember to add this export line to your ~/.bashrc or the submission script, so that EChO package is accessible by Python when you run the job.

You need to use the absolute path (you can check it by running pwd in Bash shell) for this purpose.

After these, run the following line to test:
```python
python -c 'import echo'
```
If no error occurs, EChO should have been imported into your path!


#### ASE-related settings
In the jobs submission script, besides the `PYTHONPATH`, you should also set some [environment variables](https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#environment-variables) for the `ASE` interface to `VASP`:
```bash
export PYTHONPATH=$PYTHONPATH:/xxx/echo-main
export VASPHOME=/xxx/vasp.x.x.x/bin
export VASP_COMMAND="mpirun -n $NSLOTS $VASPHOME/vasp_std"
export VASP_PP_PATH=/xxx/potentials/
```

The `VASP_COMMAND` would vary from machine to machine, depending on job scheduler and number of cores to use. You may need to adjust it yourself (and the script header, of course).


## Tutorials

See the `examples` folder.
