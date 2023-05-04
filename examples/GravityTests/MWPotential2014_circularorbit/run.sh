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
echo "Clearing existing figures."
if [ -f "circular_orbits.png" ];
then
    rm circular_orbits.png
fi

if [ -f "deviation.png" ];
then
    rm deviation.png
fi

if [ -f "deviation_from_original_data.png" ];
then
    rm deviation_from_original_data.png
fi

#Generates the initial conditions
echo "Generate initial conditions for circular orbits"
if command -v python3 &>/dev/null; then
    python3 makeIC.py
else
    python3 makeIC.py
fi

#Runs the simulation
# self gravity G, external potential g, hydro s, threads t and high verbosity v
echo "Starting the simulation... Look at output.log for the simulation details."
../../../swift --external-gravity --threads=8 MWPotential2014_circular_orbits.yml 2>&1 > output.log
echo "Simulation ended."

#Saves the plots
echo "Save plots of the circular orbits and of the errors"
if command -v python3 &>/dev/null; then
    python3 makePlots.py
else 
    python3 makePlots.py
fi
