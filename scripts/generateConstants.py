#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 16 14:27:24 2019

@author: maksbh
"""




import numpy as np
from termcolor import colored, cprint 

"""
Function to find the nearest index

Parameters
----------

array : Numpy array
value : value of energy

Returns
-------
idx : Integer
      index location corresponding to the closest location
"""

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

"""
Function to get the interpolated value

Parameters
----------

array : Numpy array
value : value of energy
nearest_id : id corresponding to the nearest value

Returns
-------
valArray : Numpy array
           array of the interpolated values
"""
def get_interpolated_value(array,value,nearest_id,energy_id):
    valArray = np.zeros(array.shape[1]);
    if(array[nearest_id][energy_id] > value):
        xp = [array[nearest_id - 1][energy_id], array[nearest_id ][energy_id]];
        for i in range(0,array.shape[1]):
            yp = [array[nearest_id - 1][i], array[nearest_id][i]];
            valArray[i] = np.interp(value,xp,yp);

    elif (array[nearest_id][energy_id] < value):
        xp = [array[nearest_id][energy_id], array[nearest_id + 1][energy_id]];
        for i in range(0,array.shape[1]):
            yp = [array[nearest_id][i], array[nearest_id + 1][i]];
            valArray[i] = np.interp(value,xp,yp);

    else:
        for i in range(0,len(valArray)):
            valArray[i] = array[nearest_id][i];

    

    return valArray;


def removeDuplicates(Data,energy_id):
    listIn = Data.tolist()
    listOut=[]
    listOut.append(listIn[0])
    currEnergy = listIn[0][energy_id]
    duplicateFound = False;
    for i in range(1,len(listIn)):
        if(listIn[i][0] == currEnergy):
            duplicateFound = True;
            continue
        else:
            listOut.append(listIn[i])
            currEnergy = listIn[i][energy_id]

    if(duplicateFound):
        cprint('Duplicates in Energy found. Removing it', 'yellow') 
    return(np.array(listOut))
    

def dump_dataVacuum(index,energy,f):
     Header = "EnergyData" + str(index) +":\n{\n";
     f.write(Header);
     Energy = "Energy = " + str(energy) + ";\n";
     f.write(Energy);
     BetaPara = "BetaPara = " + str(0.0) + ";\n";
     f.write(BetaPara);
     BetaPerp = "BetaPerp = " + str(0.0) + ";\n";
     f.write(BetaPerp);
     DeltaPara = "DeltaPara = " + str(0.0) + ";\n";
     f.write(DeltaPara);
     DeltaPerp = "DeltaPerp = " + str(0.0) + ";\n";
     f.write(DeltaPerp);
     f.write("}\n");


def dump_data(valArray,index,labelEnergy,f):
     Header = "EnergyData" + str(index) +":\n{\n";
     f.write(Header);
     Energy = "Energy = " + str(valArray[labelEnergy["Energy"]]) + ";\n";
     f.write(Energy);
     BetaPara = "BetaPara = " + str(valArray[labelEnergy["BetaPara"]]) + ";\n";
     f.write(BetaPara);
     BetaPerp = "BetaPerp = " + str(valArray[labelEnergy["BetaPerp"]]) + ";\n";
     f.write(BetaPerp);
     DeltaPara = "DeltaPara = " + str(valArray[labelEnergy["DeltaPara"]]) + ";\n";
     f.write(DeltaPara);
     DeltaPerp = "DeltaPerp = " + str(valArray[labelEnergy["DeltaPerp"]]) + ";\n";
     f.write(DeltaPerp);
     f.write("}\n");




def main(startEnergy, endEnergy,increment,dict,labelEnergy,numMaterial):
    NumEnergy = int(np.round((endEnergy - startEnergy)/increment + 1));

    for numMat in range(0,numMaterial):
        f = open("Material" + str(numMat) + ".txt", "w")
        fname = dict["Material" + str(numMat)]
        if(fname != 'vacuum'):
            Data = np.loadtxt(fname,skiprows=1);
            Data = Data[Data[:,labelEnergy["Energy"]].argsort()]
            Data = removeDuplicates(Data,labelEnergy["Energy"])
            for i in range(0,NumEnergy):
                currentEnergy = startEnergy + i* increment;
                nearest_id = find_nearest(Data[:,labelEnergy["Energy"]],currentEnergy)
                ValArray = get_interpolated_value(Data,currentEnergy,nearest_id,labelEnergy["Energy"])
                dump_data(ValArray,i,labelEnergy,f)

        else:
            for i in range(0,NumEnergy):
                energy = startEnergy + increment*i
                dump_dataVacuum(i,energy,f)
        f.close()



if __name__ == "__main__":
    startEnergy = 280.0; #Start energy
    endEnergy = 290.0;   #End  energy
    incrementEnergy = 0.1; #Increment  energy
    startAngle = 0.0; #start angle
    endAngle = 180.0; #end angle
    incrementAngle = 2.0; #increment in each angle
    numThreads = 4; #number of threads for execution
    numX = 2048; # number of voxels in X direction
    numY = 2048;# number of voxels in Y direction
    numZ = 1;# number of voxels in Z direction
    physSize = 5.0;
    rotMask = False;
    EwaldsInterpolation = 1; # 1 : Linear Interpolation 0: Nearest Neighbour
    writeVTI = False;
    windowingType = 0;

    #Files corresponding to Each material. For vacuum pass null
    dict={'Material0':'constants.txt'}


    # Label of energy to look for
    labelEnergy={"BetaPara":1,
                 "BetaPerp":2,
                 "DeltaPara":3,
                 "DeltaPerp":4,
                 "Energy":0}

#### Do not change below this
    f = open("config.txt", "w")
    f.write("StartEnergy = " + str(startEnergy) + ";\n");
    f.write("EndEnergy = " + str(endEnergy) + ";\n");
    f.write("IncrementEnergy = " + str(incrementEnergy) + ";\n");
    f.write("StartAngle = " + str(startAngle) + ";\n");
    f.write("EndAngle = " + str(endAngle) + ";\n");
    f.write("IncrementAngle = " + str(incrementAngle) + ";\n");
    f.write("NumThreads = " + str(numThreads) + ";\n");
    f.write("NumX = " + str(numX) + ";\n");
    f.write("NumY = " + str(numY) + ";\n");
    f.write("NumZ = " + str(numZ) + ";\n");
    f.write("PhysSize = " + str(physSize) + ";\n");
    f.write("RotMask = " + str(rotMask) + ";\n");
    f.write("EwaldsInterpolation= " + str(EwaldsInterpolation) + ";\n");
    f.write("WriteVTI = " + str(writeVTI) + ";\n");
    f.write("WindowingType = " + str(windowingType) + ";\n");
    f.close();



    main(startEnergy,endEnergy,incrementEnergy,dict,labelEnergy,len(dict));
