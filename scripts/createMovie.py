#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 22 16:13:54 2020

@author: maksbh
"""

import os
path = os.path.dirname(os.path.realpath(__file__))
folder = "VTI"

X = []
for file in os.listdir(folder):
    if file.endswith(".vti"):
        fname = folder + "/" + file
        X.append(fname)


#
f = open("movie.pvd", "w");
f.write("<?xml version=\"1.0\"?>\n")
f.write("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
f.write("<Collection>\n");

for i in range (0,len(X)):
    f.write("<DataSet timestep=\"" + str(i) + "\" group=\"\" part=\"0\" file=\"" + X[i] + "\"/>\n")

f.write("</Collection>\n")
f.write("</VTKFile>\n")
f.close()
