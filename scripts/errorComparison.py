#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 31 18:09:26 2019

@author: maksbh
"""

import pandas as pd
dirActual = 'sdada'
dirReference = 'dsadad'

for i in range(0,100):
    fnameActual = dirActual + '/CSV/projectionAverage' + str(i) + '.csv'
    fnameReference = dirReference + '/CSV/projectionAverage' + str(i) + '.csv'
    dfActual = pd.read_csv(fnameActual)
    dfReference = pd.read_csv(fnameReference)
    error = dfActual["projection"] - dfReference['projection']
    sortedError = error.sort_values(axis=0)
    print('Maximum Differnece = ', sortedError.iloc[-1])
    print('Minimum Differnece = ', sortedError.iloc[0])
    

    