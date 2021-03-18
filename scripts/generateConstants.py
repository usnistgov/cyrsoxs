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


def get_interpolated_value(array, value, nearest_id, energy_id):
    valArray = np.zeros(array.shape[1])
    if (array[nearest_id][energy_id] > value):
        xp = [array[nearest_id - 1][energy_id], array[nearest_id][energy_id]]
        for i in range(0, array.shape[1]):
            yp = [array[nearest_id - 1][i], array[nearest_id][i]]
            valArray[i] = np.interp(value, xp, yp)

    elif (array[nearest_id][energy_id] < value):
        xp = [array[nearest_id][energy_id], array[nearest_id + 1][energy_id]]
        for i in range(0, array.shape[1]):
            yp = [array[nearest_id][i], array[nearest_id + 1][i]]
            valArray[i] = np.interp(value, xp, yp)

    else:
        for i in range(0, len(valArray)):
            valArray[i] = array[nearest_id][i]

    return valArray


def removeDuplicates(Data, energy_id):
    listIn = Data.tolist()
    listOut = []
    listOut.append(listIn[0])
    currEnergy = listIn[0][energy_id]
    duplicateFound = False
    for i in range(1, len(listIn)):
        if (listIn[i][energy_id] == currEnergy):
            duplicateFound = True
            continue
        else:
            listOut.append(listIn[i])
            currEnergy = listIn[i][energy_id]

    if (duplicateFound):
        cprint('Duplicates in Energy found. Removing it', 'yellow')
    return (np.array(listOut))


def dump_dataVacuum(index, energy, f):
    Header = "EnergyData" + str(index) + ":\n{\n"
    f.write(Header)
    Energy = "Energy = " + str(energy) + ";\n"
    f.write(Energy)
    BetaPara = "BetaPara = " + str(0.0) + ";\n"
    f.write(BetaPara)
    BetaPerp = "BetaPerp = " + str(0.0) + ";\n"
    f.write(BetaPerp)
    DeltaPara = "DeltaPara = " + str(0.0) + ";\n"
    f.write(DeltaPara)
    DeltaPerp = "DeltaPerp = " + str(0.0) + ";\n"
    f.write(DeltaPerp)
    f.write("}\n")


def dump_data(valArray, index, labelEnergy, f):
    Header = "EnergyData" + str(index) + ":\n{\n";
    f.write(Header)
    Energy = "Energy = " + str(valArray[labelEnergy["Energy"]]) + ";\n"
    f.write(Energy)
    BetaPara = "BetaPara = " + str(valArray[labelEnergy["BetaPara"]]) + ";\n"
    f.write(BetaPara)
    BetaPerp = "BetaPerp = " + str(valArray[labelEnergy["BetaPerp"]]) + ";\n"
    f.write(BetaPerp)
    DeltaPara = "DeltaPara = " + str(valArray[labelEnergy["DeltaPara"]]) + ";\n"
    f.write(DeltaPara)
    DeltaPerp = "DeltaPerp = " + str(valArray[labelEnergy["DeltaPerp"]]) + ";\n"
    f.write(DeltaPerp)
    f.write("}\n")


def writeList(name: str, value: list, file):
    valStr: str = name + "["
    for i in range(len(value) - 1):
        valStr = valStr + str(value[i]) + ","
    valStr = valStr + str(value[len(value) - 1])
    file.write(valStr + "];\n")


def main(energies, dict, labelEnergy, numMaterial):
    NumEnergy = len(energies)

    for numMat in range(0, numMaterial):
        f = open("Material" + str(numMat) + ".txt", "w")
        fname = dict["Material" + str(numMat)]
        if (fname != 'vacuum'):
            Data = np.loadtxt(fname, skiprows=1)
            Data = Data[Data[:, labelEnergy["Energy"]].argsort()]
            Data = removeDuplicates(Data, labelEnergy["Energy"])
            for i in range(0, NumEnergy):
                currentEnergy = energies[i]
                nearest_id = find_nearest(Data[:, labelEnergy["Energy"]], currentEnergy)
                ValArray = get_interpolated_value(Data, currentEnergy, nearest_id, labelEnergy["Energy"])
                dump_data(ValArray, i, labelEnergy, f)

        else:
            for i in range(0, NumEnergy):
                energy = startEnergy + increment * i
                dump_dataVacuum(i, energy, f)
        f.close()


if __name__ == "__main__":
    energies: list = [280.0, 285.0, 281.0]
    eAngleRotation: list = [0.0, 2.0, 180.0]  # [start : increment: end]
    numThreads = 4  # number of threads for execution
    numX = 2048  # number of voxels in X direction
    numY = 2048  # number of voxels in Y direction
    numZ = 1  # number of voxels in Z direction
    kRotationType = 0 # 0: No rotation 1: kRotation
    kAngleRotation: list = [0.0, 2.0, 180.0]  # [start : increment: end]
    physSize = 5.0
    rotMask = False
    EwaldsInterpolation = 1  # 1 : Linear Interpolation 0: Nearest Neighbour
    writeVTI = False
    windowingType = 0
    morphologyType = 0  # 0: Euler angles 1: Vector Morphology 2: Spherical coordinate
    scatterApproach = 0  # 0 : Partial (Default) 1: Full

    # Files corresponding to Each material. For vacuum pass vacuum
    dict = {'Material0': '../OpticalConstants/PEOlig2018.txt'}

    # Label of energy to look for
    labelEnergy = {"BetaPara": 0,
                   "BetaPerp": 1,
                   "DeltaPara": 2,
                   "DeltaPerp": 3,
                   "Energy": 6}

    #### Do not change below this
    f = open("config.txt", "w")
    writeList("Energies=", value=energies, file=f)
    writeList("EAngleRotation=", value=eAngleRotation, file=f)
    f.write("NumThreads = " + str(numThreads) + ";\n")
    f.write("NumX = " + str(numX) + ";\n")
    f.write("NumY = " + str(numY) + ";\n")
    f.write("NumZ = " + str(numZ) + ";\n")
    f.write("PhysSize = " + str(physSize) + ";\n")
    f.write("RotMask = " + str(rotMask) + ";\n")
    f.write("EwaldsInterpolation= " + str(EwaldsInterpolation) + ";\n")
    f.write("WriteVTI = " + str(writeVTI) + ";\n")
    f.write("WindowingType = " + str(windowingType) + ";\n")
    f.write("MorphologyType = " + str(morphologyType) + ";\n")
    f.write("ScatterApproach = " + str(scatterApproach) + ";\n")
    f.write("KRotationType = " + str(kRotationType) + ";\n")
    writeList("KAngleRotation=", value=kAngleRotation, file=f)
    f.close()

    main(energies, dict, labelEnergy, len(dict))
