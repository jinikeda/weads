#!/bin/bash
python preprocessing.py ;
python hydrodynamics.py ;
python interpdatums.py ;
python ecology.py ;
python postprocessing.py ;
wait

