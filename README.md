Spinner is a magnetochemical software package.
It uses the Spin-Hamiltonian approach and statistical physics to describe the magnetism of molecular systems.
At the moment, Spinner works only with the temperature dependence of the
magnetic susceptibility.

## Build

### 1. Clone Spinner:

```shell
git clone https://github.com/ruthenium96/spinner.git
cd spinner
```

### 2. Install dependencies

#### Manual installation

Install following dependencies:

- `Armadillo` or `Eigen3`;
- `OpenMP`;
- Fortran compiler (required for wignerSymbols);
- `GoogleTest` (required for tests, optional).

#### via Docker:

Install [Docker](https://www.docker.com/) and run:

```shell
chmod +x docker.sh
sudo ./docker.sh init
sudo ./docker.sh shell # or sudo ./docker.sh shell_restart
cd project
```

### 3. Build the project:

```shell
mkdir build
cd build
cmake ..
make # or cmake --build .
```
## Usage

```shell
spinner_main file.yml
```
where `file.yml` is the input file containing all the information needed for the calculation.

The input file `file.yml` currently consists of four mandatory blocks. 

### `model_input`

To define model(s) it is necessary to specify a set of multiplicities and a set of spin Hamiltonian parameters. There are three modes of operation of the program: `single`, as well as `trajectory` and `scan`, the last two potentially useful for the construction of residual surfaces.

One can specify `centers` using either a list or a key `all`.

```yml
model_input:
  multiplicities: [4, 4, 4, 4, 4, 4, 4, 4]
  mode: single
  parameters:
    - name: J
      type: J
      value: +100.0
      fixed: false
      pairs: [[0, 1], [1, 2], [2, 3], [4, 5], [5, 6], [6, 7],
              [0, 4], [1, 5], [2, 6], [3, 7]]
    - name: g
      type: g_factor
      value: 2.0
      fixed: true
      centers: all
    - name: theta
      type: theta
      value: -5.0
      fixed: false
```

### `optimizations`

There are three different `mode`s of `optimization`: `none`, `custom` and `auto` (not implemented yet).

```yml
optimizations:
  mode: none
```

In the `custom` mode, it is nessesary to specify initial basis, either `lex` or `ito`. it is also possible to specify different optimizations: `tz_sorter`, `tsquared_sorter`, `positive_tz_eliminator`, `nonminimal_tz_eliminator`, `s2_transformer` and `symmetrizer`. `symmetrizer` should be a list of commuting groups, to specify group one need to specify name of group (`S2` or `Dihedral` + its order) and list of generators of this group.

```yml
optimizations:
  mode: custom
  custom:
    basis: ito
    tz_sorter:
    tsquared_sorter:
    nonminimal_tz_eliminator:
    symmetrizer:
      - group_name: S2
        generators: [[6, 7, 8, 3, 4, 5, 0, 1, 2]]
      - group_name: S2
        generators: [[2, 1, 0, 5, 4, 3, 8, 7, 6]]
```

### `job`
This block controls the type of work Spinner has to do. The most important key is `mode`, which can have two values: `simulation` and `fit`.

#### `simulation`
In this case, it is necessary to specify only those temperatures for which it is necessary to construct a theoretical temperature dependence. This can be done either by using a list of required temperatures:
```yml
job:
  mode: simulation
  simulation:
    temperatures: [1, 2, 5, 100, 150, 200]
```
or by specifying the start and end temperatures and the step between the temperatures:
```yml
job:
  mode: simulation
  simulation:
    temperatures: !range [1.0, 300.0, 2]
```

#### `fit`

In this case, it is necessary to determine the nonlinear regression method (`optim_nm` or `stl_bfgs`) and to specify the experimental data for which the inverse problem must be solved. To specify the experimental data, it is necessary to specify 
- the ratio of the number of centres in the experiment to the number of centres in the model, 
- the dimensionality of the experimental data (`mu_in_bohr_magnetons`,`mu_squared_in_bohr_magnetons_squared`, `chiT_in_cm_cubed_kelvin_per_mol`, `chi_in_cm_cubed_per_mol`), 
- the method of weighting the experimental data (`per_point`, `per_interval`),
- the data themselves, either in the form of a list of pairs (temperature, value) or as a path to a file with data in two columns.

```yml
job:
  mode: fit
  fit:
    solver: optim_nm
    experiment:
      data: !file experiment.exp
      dimension: mu_in_bohr_magnetons
      ratio: 2.0
      weights: per_interval
```
```yml
job:
  mode: fit
  fit:
    solver: optim_nm
    experiment:
      data: [[1, 1], [100, 0.5], [200, 0.25]]
      dimension: mu_in_bohr_magnetons
      ratio: 2.0
      weights: per_point
```

### `control`

This block consists of three mandatory keys: the level of print detail (`trace`, `debug`, `verbose`, `detailed`, `basic`, `error`, and `off`), the precision of the floating-point numbers (`single` or `double`), and the linear algebra package to be used for the calculation (`arma` or `eigen`).
Example:
```yml
control:
  print_level: trace
  dense_precision: double
  dense_algebra_package: arma
```