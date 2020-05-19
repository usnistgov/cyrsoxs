GPU enabled RSoXS simulation (Release 1.0 - Beta version)
====================================
## Generating config and constant files

In order to run RSoXs, first we need to generate the confiig file.
To generate the config file, you need to run the python script `generateConstants.py` 
from the scripts folder.

Here are the description of the individual parameters that you would need to change
in `__main__`.
```
startEnergy = 280.0; #Start energy
endEnergy = 290.0;   #End  energy
incrementEnergy = 0.1; #Increment in  energy
startAngle = 0.0; #start angle
endAngle = 180.0; #end angle
incrementAngle = 2.0; #increment in each angle
numThreads = 4; number of threads for execution
numX = 2048; # number of voxels in X direction
numY = 2048;# number of voxels in Y direction
numZ = 1;# number of voxels in Z direction
physSize = 5.0; #Physical size
RotMask = False; #Default: False
EwaldsInterpolation = 1; # 1 : Linear Interpolation (default) 0: Nearest Neighbour 
WriteVTI = False; # Default : False
WindowingType = 0; # 0: None (Default) 1: Hanning 
``` 

This code also generate the optical constants for each Energy level
by interpolating from the files provided.

```
dict={'Material0':'Filename for Material Constants 1',
      'Material1': 'Filename for Material Constants 2',
      'Material2':'Filename for Material Constants 3',
      'Material3':'Filename for Material Constants 4'}
```
As an example:
```
dict={'Material0':'../OpticalConstants/PEOlig2018.txt',
      'Material1': '../OpticalConstants/PEOlig2018.txt',
      'Material2':'../OpticalConstants/PEOlig2018.txt',
      'Material3':'vacuum'}
```

Here, Material0, Material1, Material2 constants are read from the 
file `PEOlig2018.txt`. The Material3 property is set to vacuum.

Also, you need to provide the individual columnID.

```
labelEnergy={"BetaPara":0,
             "BetaPerp":1,
             "DeltaPara":2,
             "DeltaPerp":3,
             "Energy":6}
```  

This basically  says that the `0` column of the opticalConstants file corresponds to `BetaPara` , 1 
corresponds to `BetaPerp` and so on.

The python script can be ran directly by:
```
python generateConstants.py
``` 

Once the script has successfully completed, it will generate the files `config.txt` and `Material0.txt` ,
 `Material1.txt` and so on for each individual material. 
 
 
Running RSoXS
=============

Copy all the generated files to run directory.

 In order to run Cy-RSOXS you need to execute the following command
from the run directory.

```
./$(PATH_TO_CY-RSoXS_BUILD_DIR)/Cy-RSoXS  $(PATH_TO_HDF5_FILE)
```

The output will be generated in the folder named `Projection`that will be
created in the run directory. The output are generated in `.vti` / `.hdf5` format which
can be visualized using [Paraview](https://www.paraview.org/) or [Visit](https://wci.llnl.gov/simulation/computer-codes/visit/).
