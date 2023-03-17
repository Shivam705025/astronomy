#!/bin/bash

rm output/output_*
rm xyplot.png xzplot.png yzplot.png
# if [ ! -e circular_orbits_MW.hdf5 ]
# then
echo "Generate initial conditions for circular orbits"
if command -v python3 &>/dev/null; then
    python3 makeIC.py
else
    python3 makeIC.py
fi

# fi

# self gravity G, external potential g, hydro s, threads t and high verbosity v
../../../swift --external-gravity --threads=8 nfw_mn_psc_circular_orbits.yml 2>&1 | tee output.log

echo "Save plots of the circular orbits"
if command -v python3 &>/dev/null; then
    python3 plotprog.py
else 
    python3 plotprog.py
fi
