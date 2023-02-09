# GPU enabled RSoXS simulation (1.1.6.0)

CyRSoXS prints output in `.h5` and  `.vti` format. Both these file format can
be easily visualized in Paraview. The output in `.h5` is stored in the directory
`H5` and the `.vti` files are stored in the directory name `VTI`.

## Installing Paraview

Paraview can be downloaded from [here](https://www.paraview.org/).

## Generating a movie

1. We provide a python script `createMovie.py` located in scripts folder.
 Run `python createMovie.py` outside of VTI folder. It will generate a Paraview Data file `movie.pvd`.
  
2. Once `movie.pvd` is created, open Paraview and click on `File->LoadState`. Load the State file `movie.pvsm`
located under scripts. This will prompt you to provide the location of `movie.pvd` to load.

3. Once you have loaded `movie.pvd` into Paraview, you can play the movie by clicking on the play icon.
