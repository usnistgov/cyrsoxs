GPU enabled RSoXS simulation (0.9.0)
====================================

Installing Docker
=================

On Ubuntu you can use the command to install docker:

```
sudo apt-get install docker
sudo apt-get install docker.io 
```

Launching Docker
================
`sudo docker pull maksbh/cy-rsoxs:tagname`

To get the imageID run:

`sudo docker images`

It will show the output as:

```
REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
maksbh/cy-rsoxs     latest              655455b309d5        17 minutes ago      4.77GB
```

Installing Cy-RSoXS
====================

First you need to run the docker interactively:

`sudo docker run -it $IMAGE_ID`


Then run the following command:
```
cd; ./configure;
```

It will ask for your bitbucket username and password. The code will automatically compile and you can run.

Note that you might need to change the configure file. Here are the contents of configuration file:

```
#clone the git repository
git clone https://bitbucket.org/baskargroup/cy-rsoxs.git

# Go into the directory
cd cy-rsoxs;

mkdir build; cd build;

#Cmake option.Might need to change that.
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_BUILD_TYPE=Release -DHDF5_DIR=/root/Dependencies/CMake-hdf5-1.10.5/build/_CPack_Packages/Linux/TGZ/HDF5-1.10.5-Linux/HDF_Group/HDF5/1.10.5/share/cmake/hdf5 -DDLEVEL2=Yes -DNUM_MATERIAL=4

#Make the executable
make
```



Contributors
============
* Kumar Saurabh
* Adarsh Krishnamurthy
* Baskar Ganapathysubramanian
* Eliot Gann
* Dean Delongchamp
* Michael Chabinyc

Acknowledgement
===============
We thank ONR MURI Center for Self-Assembled Organic Electronics for providing the support for this work.

Contact
=======
Questions and comments should be sent to Dr. Baskar Ganapathysubramanian at [baskarg@iastate.edu](mailto:baskarg@iastate.edu) or Dr.  Adarsh Krishnamurthy [adarsh@iastate.edu](mailto:adarsh@iastate.edu).
