# Crystal PBM Models

Custom OpenQBMM population balance models and solvers for crystallization.

## Prerequisites

- OpenFOAM v2412. USed this because OpenQBMM webpage suggests this is the last one supported. However, later I found out that web-page is outdated, and in their git they update regularly. Latest is v2606. 
- OpenQBMM (built from source)

## Installation

### 1. Install OpenFOAM

Install OpenFOAM v2412 following the official instructions.

### 2. Build OpenQBMM

Clone and build OpenQBMM:

```bash
git clone <OpenQBMM repository>
cd OpenQBMM
./Allwmake
```

### 3. Configure the environment

Add the location of the OpenQBMM source tree to your `~/.bashrc`, e.g.:

```bash
export QBMM_DIR=$HOME/OpenFOAM/username-2412/OpenQBMM
```

Reload your shell:

```bash
source ~/.bashrc
```

### 4. Build this repository

```bash
./Allwmake
```

## Using solvers

Check out examples folder!
