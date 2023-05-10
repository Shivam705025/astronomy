#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
# parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("output")
argparser.add_argument("--ncell", "-n", type=int, default=500)
args = argparser.parse_args()

the_statistics=np.transpose(np.loadtxt("statistics.txt"))

Time = the_statistics[1]
E_kin = the_statistics[13]
E_int = the_statistics[14]
E_mag = the_statistics[34]
E_tot=E_kin+E_int+E_mag

fig, ax = plt.subplots(1, 2, sharex=True, figsize=(10, 5))
ax[0].plot(Time, E_kin/E_tot,label="E_kin")
ax[0].plot(Time, E_int/E_tot,label="E_int")
ax[0].plot(Time, E_mag/E_tot,label="E_mag")
ax[0].plot(Time, E_tot/E_tot,label="E_tot",color="Black",linestyle="dashed")
ax[0].set_xlabel("Time [s]")
ax[0].set_ylabel("E/ E_tot")
ax[0].legend(loc="best")

ax[1].plot(Time, E_tot/E_tot[0],label="E(t)")
ax[1].plot(Time, E_tot/E_tot,label="Energy conservation",color="Black",linestyle="dashed")
ax[1].set_xlabel("Time [s]")
ax[1].set_ylabel("E_tot(t) / E_tot(0)")
ax[1].legend(loc="best")
plt.tight_layout()
plt.savefig(args.output+'/'+'E_change.png', dpi=100)
