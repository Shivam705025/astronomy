#!/bin/bash

#Creates a directory for the outputs
DIR=output
if [ -d "$DIR" ];
then
    echo "$DIR directory exists. Its content will be removed."
    rm output/output_*
else
    echo "$DIR directory does not exists. It will be created."
    mkdir $DIR
fi

#Clears the previous figures
rm xyplot.png xzplot.png yzplot.png

#Generates the initial conditions
echo "Generate initial conditions for circular orbits"
if command -v python3 &>/dev/null; then
    python3 makeIC.py
else
    python3 makeIC.py
fi

#Runs the simulation
# self gravity G, external potential g, hydro s, threads t and high verbosity v
../../../swift --external-gravity --threads=8 nfw_mn_psc_circular_orbits.yml 2>&1 | tee output.log

#Saves the plots
echo "Save plots of the circular orbits"
if command -v python3 &>/dev/null; then
    python3 plotprog.py
else 
    python3 plotprog.py
fi
