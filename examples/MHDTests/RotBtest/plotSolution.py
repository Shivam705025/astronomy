#!/usr/bin/env python3

"""
plotSolution.py

Plot the density, pressure, divB, velocity components and magnetic field
components for the Brio Wu snapshot input and create a figure with the
given output name.

Usage:
  python3 plotSolution.py SNAPSHOT OUTPUT [--ncell NCELL]
where NCELL is the number of cells to use for the HLL Riemann solver reference
solution (default: 1000).

Also plots the "exact" solution that is hard-coded in exact_solution.py.
"""

import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import argparse
# parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("input")
argparser.add_argument("output")
argparser.add_argument("--ncell", "-n", type=int, default=500)
args = argparser.parse_args()

# read the variables of interest from the snapshot file
gamma = None
boxsize = None
t = None
x = None
rho = None
v = None
P = None
B = None
divB = None
rotB = None

with h5py.File(args.input, "r") as handle:
    gamma = handle["HydroScheme"].attrs["Adiabatic index"][0]
    boxsize = handle["Header"].attrs["BoxSize"][0]
    #print(handle["/PartType0"].keys())
    #print(handle["/HydroScheme"].attrs.keys())
    t = handle["Header"].attrs["Time"][0]

    x = handle["PartType0/Coordinates"][:, 0]
    y = handle["PartType0/Coordinates"][:, 1]
    z = handle["PartType0/Coordinates"][:, 2]

    rho = handle["PartType0/Densities"][:]
    h = handle["PartType0/SmoothingLengths"][:]
    v = handle["PartType0/Velocities"][:]
    P = handle["PartType0/Pressures"][:]
    B = handle["PartType0/MagneticFluxDensity"][:]
    divB = handle["PartType0/MagneticDivergence"][:]
    monB = handle["PartType0/MonopolePartB"][:]
    #divBest = handle["PartType0/BmonRatio"][:]
    git = handle["Code"].attrs["Git Revision"]
    gitBranch = handle["Code"].attrs["Git Branch"]
    scheme = handle["/HydroScheme"].attrs["Scheme"]
    kernel = handle["/HydroScheme"].attrs["Kernel function"]
    neighbours = handle["/HydroScheme"].attrs["Kernel target N_ngb"]
    mu0 = handle["/PhysicalConstants/InternalUnits"].attrs["vacuum_permeability"]
   # dedhyp = handle["/HydroScheme"].attrs["Dedner Hyperbolic Constant"]
   # dedpar = handle["/HydroScheme"].attrs["Dedner Parabolic Constant"]
   # mhdeta = handle["/HydroScheme"].attrs["Diffusion Eta"]
   # mhdflavour = handle["/HydroScheme"].attrs["MHD Flavour"]
bb = np.sqrt(B[:,0]*B[:,0]+B[:,1]*B[:,1]+B[:,2]*B[:,2])
x -= boxsize * 0.5
#Abs_F_L=np.sqrt((rotB[:,1]*B[:,2]-rotB[:,2]*B[:,1])**2+(rotB[:,2]*B[:,0]-rotB[:,0]*B[:,2])**2+(rotB[:,0]*B[:,1]-rotB[:,1]*B[:,0])**2)

#Cut_edge
d_edg=0.1
maskCube=np.where((x>np.min(x)+d_edg) & (x<np.max(x)-d_edg) & (y>np.min(y)+d_edg) & (y<np.max(y)-d_edg) & (z>np.min(z)+d_edg)& (z<np.max(z)-d_edg))
print(f'{len(x)} and {len(x[maskCube])}')
x=x[maskCube]
y=y[maskCube]
z=z[maskCube]
rho=rho[maskCube]
h=h[maskCube]
v=v[maskCube]
P=P[maskCube]
B=B[maskCube]
divB=divB[maskCube]
monB=monB[maskCube]
#divBest=divBest[maskCube]

print(np.max(np.abs(monB[:,0])),len(monB[:,0]))
# plot everythinplots(3, 3, sharex=True, figsize=(10, 9))
# plot everything
fig, ax = pl.subplots(3, 3, sharex=True, figsize=(10, 9))
mms=1
aalpha=0.3

ax[0][0].plot(x, rho, ".",markersize=mms,alpha=aalpha)
ax[0][1].plot(x, P, ".",markersize=mms,alpha=aalpha)
ax[0][2].plot(x, divB, ".",markersize=mms,alpha=aalpha,label="divB")
monB=monB*100
#ax[0][2].plot(x,divBest,".",markersize=mms,alpha=aalpha,label="|divBest|")
for i in range(3):
    ax[1][i].plot(x, v[:, i], ".",markersize=mms,alpha=aalpha)
    
ax[2][0].plot(x, B[:, 0], ".", label="SWIFT",markersize=mms,alpha=aalpha)
ax[2][0].plot(x,0.013*monB[:,0],".",markersize=mms,alpha=aalpha,label="|B_mon_x|")
ax[2][1].plot(x, B[:, 1], ".",markersize=mms,alpha=aalpha)
ax[0][0].set_ylabel("rho")
ax[0][1].set_ylabel("P")
ax[0][2].set_ylabel("Quantity")
ax[1][0].set_ylabel("vx")
ax[1][1].set_ylabel("vy")
ax[1][2].set_ylabel("vz")
ax[2][0].set_ylabel("Bx")
ax[2][1].set_ylabel("By")
#ax[2][2].set_ylabel("Bz")
#ax[0][2].set_yscale("log")
ax[2][0].legend(loc="best")
ax[0][2].legend(loc="best")

# Information -------------------------------------
#plt.subplot(236, frameon=False)

text_fontsize = 7

ax[2][2].text(-0.45, 0.8,
    "BrioWu shock Tube with $\\gamma=%.3f$ in 3D\nat $t=%.2f$" % (gamma, t),
    fontsize=text_fontsize,
)
ax[2][2].plot([-0.45, 0.45], [0.62, 0.62], "k-", lw=1)
ax[2][2].text(-0.45, 0.5, "$SWIFT$ %s" % git.decode("utf-8"), fontsize=text_fontsize)
ax[2][2].text(-0.45, 0.4, "$Branch$ %s" % gitBranch.decode("utf-8"), fontsize=text_fontsize)
ax[2][2].text(-0.45, 0.3, scheme.decode("utf-8"), fontsize=text_fontsize)
ax[2][2].text(-0.45, 0.2, kernel.decode("utf-8"), fontsize=text_fontsize)
ax[2][2].text(-0.45, 0.1, "$%.2f$ neighbours" % (neighbours),fontsize=text_fontsize)
ax[2][2].plot([-0.45, 0.45], [0.0, 0.0], "k-", lw=1)
ax[2][2].text(-0.45, -0.1, "$\\mu_0:%.2f$ " % (mu0),fontsize=text_fontsize)
#ax[2][2].text(-0.45, -0.2, "$Difussion_\\eta:%.2f$ " % (mhdeta),fontsize=text_fontsize)
#ax[2][2].text(-0.45, -0.3, "$Dedner: Hyp=%.2f // Par=%.2f$" % (dedhyp,dedpar),fontsize=text_fontsize)
ax[2][2].plot([-0.45, 0.45], [-0.4, -0.4], "k-", lw=1)
#ax[2][2].text(-0.45, -0.5, "$Flavour: $ %s" % mhdflavour.decode("utf-8")[0:30], fontsize=text_fontsize)

for a in ax[2]:
    a.set_xlabel("x")

# make sure all velocity and magnetic fields plots use the same y axis
# this puts the noise on constant values into context
#for a in [ax[1][0], ax[1][2]]:
#    a.set_ylim(*ax[1][1].get_ylim())
#for a in [ax[2][0], ax[2][2]]:
#    a.set_ylim(*ax[2][1].get_ylim())
ax[0][0].set_ylim(0.0,1.1)
ax[0][1].set_ylim(0.0,1.1)
#ax[0][2].set_ylim(-0.2,0.2)
ax[1][0].set_ylim(-2.0,1.0)
ax[1][1].set_ylim(-2.0,1.0)
ax[1][2].set_ylim(-2.0,1.0)
ax[2][0].set_ylim(-1.1,1.1)
ax[2][1].set_ylim(-1.1,1.1)
ax[2][2].set_ylim(-1.1,1.1)
#ax[0][2].set_ylim(0.5,2)

ax[2][2].tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
ax[2][2].set_xlabel("")
ax[2][2].plot(frameon=False)

# only plot the relevant part of the solution
ax[0][0].set_xlim(-0.5, 0.5)

# add the time as a title 
ax[0][1].set_title(f"t={t:.2e}")

# save the figure
pl.tight_layout()
pl.savefig(args.output, dpi=300)
