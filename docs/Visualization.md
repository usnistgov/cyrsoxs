GPU enabled RSoXS simulation (1.0.0 - Beta)
====================================

Cy-RSoXS prints output in `.h5` and  `.vti` format. Both these file format can 
be easily visualized in Paraview. The output in `.h5` is stored in the directory 
`H5` and the `.vti` files are stored in the directory name `VTI`.

Installing Paraview
===================
The Paraview can be downloaded from [here](https://www.paraview.org/).


Generating movie
================
1. We provide a python script `createMovie.py` located in scripts folder.
 Run `python createMovie.py` outside of VTI folder. It will generate a file `movie.pvd`.
  
2. Once, `movie.pvd` is created. Open Paraview and click on `File->LoadState`. Load `movie.pvsm`
located under scripts. This will ask to give the location of `movie.pvd`. 

3. Then you can play the movie by clicking on the play icon.   