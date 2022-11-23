# ISI
Code and sample data for analysis of PET scans with in-scan drug intervention. 

This repository contains code that fits a PET compartment model with time-varying occupancy to a PET displacement scan. Two different solutions to the model are available, one Euler Forward-based numerical solution, and one analytical solution where the the change in occupancy is modelled to happen instantaneously in a single step. 

The function 'isi.m' acts as a wrapper, and calls on either 'isiFit_numerical.m'  or 'isiFit_singlestep.m', depending on what solution is specified. 'isi.m' has two inputs: the structure 'Data', holding information about the PET experiment and arterial input function, and the structure 'Settings', with options on how to run the code. See the 'isi.m' file for information on how to use it. Currently, the code is only implemented for 1TC kinetics, and tested on 11C-UCB-J data from humans and pigs.

The code has been tested on Matlab 2021a (version 9.10)

The file 'isiTestData.mat' contains simulated test data, stored in the structure 'Data', as well as a structure 'Settings', which can be used to test the code. 
