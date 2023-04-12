#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse 

# parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("folder")
args = argparser.parse_args()


folder=args.folder
file_to_load='/statistics.txt'
out_image='/E(t).png'
the_statistics=np.transpose(np.loadtxt(folder+file_to_load))

Time = the_statistics[1]
E_kin = the_statistics[13]
E_int = the_statistics[14]
E_mag = the_statistics[34]
E_tot=E_kin+E_int+E_mag

divB_err=the_statistics[35]

n_of_plots_to_display=3
fig_size=5
fig, ax = plt.subplots(1, n_of_plots_to_display, sharex=True, figsize=(n_of_plots_to_display*fig_size, fig_size))
ax[0].plot(Time, E_kin/E_tot,label="E_kin")
ax[0].plot(Time, E_int/E_tot,label="E_int")
ax[0].plot(Time, E_mag/E_tot,label="E_mag")
ax[0].plot(Time, E_tot/E_tot,label="E_tot",color="Black",linestyle="dashed")
ax[0].set_xlabel("Time [s]")
ax[0].set_ylabel("E/ E_tot")
ax[0].legend(loc="best")

ax[1].plot(Time, E_tot/E_tot[0]-1,label="E(t)")
ax[1].plot(Time, 0*E_tot,label="Energy conservation",color="Black",linestyle="dashed")
ax[1].set_xlabel("Time [s]")
ax[1].set_ylabel("E_tot(t) / E_tot(0)-1")
ax[1].legend(loc="best")

ax[2].plot(Time, divB_err, label = "divB error")
ax[2].plot(Time,0*divB_err, label = "zero error", linestyle="dashed")
ax[2].set_xlabel("Time [s]")
ax[2].set_ylabel("div B err")
ax[2].legend(loc="best")

plt.tight_layout()
plt.savefig(folder+out_image, dpi=100)
