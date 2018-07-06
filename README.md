# WRF-rLake
Code for the WRF-rLake model and data necessary for running it.

A quick glimpse of the repository:
-> WRF-rLake.f90 is the source code of the WRF-rLake model
-> run.py is the main code the run WRF-rLake under python
-> data contains several txt documents for initialization, foring the model, and for evaluation.

#################
1. the WRF-rLake model is writen in Fortran 90, but it can be run by Phython through hybrid programming
2. to run WRF-rLake under Python, please do as follows (on a Mac):
   -> download WRF-rLake.f90, run.py, and the data folder to one same routine
   -> open terminal, and make lake.so by typing in: f2py -c -m lake WRF-rLake.f90
   -> run run.py by typing in: python run.py
   -> then you will get two output documents named: output8.txt and output20.txt
