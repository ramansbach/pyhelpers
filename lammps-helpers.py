# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 09:32:24 2016

@author: rachael
In this file, we have some simple setup scripts to initialize LAMMPS data files
"""

def ljfromgro(groname,outname):
    #this function converts from a gro file containing a single atom type to a 
    #lammps data file of LJ type

