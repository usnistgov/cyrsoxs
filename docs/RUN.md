# Running CyRSoXS (1.1.5.2)

For more information on the structure of each input file, please refer to [docs/DATA.md](docs/Data.md).

## With Pybind support

Refer to the [Jupyter-notebook](../notebook/CyRSoXS.ipynb) for the instructions.

## Without Pybind support

### Generating config and constant files

In order to run CyRSoXS, first we need to generate the config file.
To generate the config file, you need to run the python script `generateConstants.py`
from the scripts folder.

Default parameters are written on lines 155-166 and can be changed if desired.

```console
#Required options
caseType = 0
energies: list = [280.0, 285.0, 281.0]
eAngleRotation: list = [0.0, 2.0, 180.0]  # [start : increment: end]
morphologyType = 0  # 0: Euler angles 1: Vector Morphology

#Optional options
numThreads = 4 # number of threads for execution
RotMask = False #Default: False
EwaldsInterpolation = 1 # 1 : Linear Interpolation (default) 0: Nearest Neighbour
WindowingType = 0 # 0: None (Default) 1: Hanning
scatterApproach = 0  # 0 : Partial (Default) 1: Full
Algorithm=1
DumpMorphology=True
MaxStreams = 1
```

This code also generates the optical constants for each Energy level
by interpolating from the files provided.

```console
dict={'Material0':'Filename for Material Constants 1',
      'Material1': 'Filename for Material Constants 2',
      'Material2':'Filename for Material Constants 3',
      'Material3':'Filename for Material Constants 4'}
```

As an example:

```console
dict={'Material0':'../OpticalConstants/PEOlig2018.txt',
      'Material1': '../OpticalConstants/PEOlig2018.txt',
      'Material2':'../OpticalConstants/PEOlig2018.txt',
      'Material3':'vacuum'}
```

Here, Material0, Material1, Material2 constants are read from the
file `PEOlig2018.txt`. The Material3 property is set to vacuum.

Also, you need to provide the individual columnID.

```console
labelEnergy={"BetaPara":0,
             "BetaPerp":1,
             "DeltaPara":2,
             "DeltaPerp":3,
             "Energy":6}
```  

This  says that the `0` column of the opticalConstants file corresponds to `BetaPara` , 1
corresponds to `BetaPerp` and so on.

The python script can be ran directly by:

```bash
python generateConstants.py
```

Once the script has successfully completed, it will generate the files `config.txt` and `Material1.txt` ,
 `Material2.txt` and so on for each individual material.

## Running CyRSoXS

Copy all the generated files to run directory.

 In order to run CyRSoXS you need to execute the following command
from the run directory.

```bash
./$(PATH_TO_CyRSoXS_BUILD_DIR)/CyRSoXS  $(PATH_TO_HDF5_FILE)
```

The output will be generated in the folder named `HDF5` for HDF5 files and `VTI` for VTI files (if dumped)
that will be created in the run directory. The output are generated in `.vti` / `.hdf5` format which
can be visualized using [Paraview](https://www.paraview.org/) or [Visit](https://wci.llnl.gov/simulation/computer-codes/visit/).
