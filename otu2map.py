# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 11:16:43 2015

@author: astride
"""

from __future__ import print_function
import pandas as pd
import sys
from os.path import expanduser

datafile = expanduser(sys.argv[1])

tabla = pd.read_table("/home/astride/MEGA/Microbiota_CX/PCA/gg_otu_table.mappingfile.csv", sep=",", index_col=0, header=0)
indices = tabla.index
muestra = tabla.columns

for x in indices:
    for y in muestra:
         if y == "#SampleID":
             continue
         
         valor = tabla.loc[x][y]
         if valor > 0:
             print(str(x) + "\t" + str(y) + "\t" + str(valor))
            