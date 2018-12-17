# WRF-rLake

Code for the WRF-rLake model and data necessary for running it.
The WRF-rLake model is writen in Fortran 90, but it can be run by Python through hybrid programming

## Structure
This repo consists of the following parts:
1. `WRF-rLake.f90`: the source code of the WRF-rLake model
2. `run.py`: the main code the run WRF-rLake under python
3. `data`: contains several txt documents for initialization, foring the model, and for evaluation.

## Usage
to run WRF-rLake under Python, please do as follows (tested working on a Mac):
1.  download `WRF-rLake.f90`, `run.py`, and the `data` directory and put them in the same folder. 

2. open terminal, and make lake.so by typing in: 
```
f2py -c -m lake WRF-rLake.f90
```
3. run `run.py`: 
```
python run.py
```

4. then you will get two output documents named: `output8.txt` and `output20.txt`
