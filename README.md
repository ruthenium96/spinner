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
