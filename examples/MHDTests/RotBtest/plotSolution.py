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
    curlB = handle["PartType0/CurlB"][:]
    mon_est_B = handle["PartType0/MonopolePartB"][:]
    sol_est_B = handle["PartType0/SolenoidalPartB"][:]
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
y -= boxsize * 0.5
z -= boxsize * 0.5
#Abs_F_L=np.sqrt((rotB[:,1]*B[:,2]-rotB[:,2]*B[:,1])**2+(rotB[:,2]*B[:,0]-rotB[:,0]*B[:,2])**2+(rotB[:,0]*B[:,1]-rotB[:,1]*B[:,0])**2)

#Cut_edge
d_edg=0.1
maskCube=np.where((x>np.min(x)+d_edg) & (x<np.max(x)-d_edg) & (y>np.min(y)+d_edg) & (y<np.max(y)-d_edg) & (z>np.min(z)+d_edg)& (z<np.max(z)-d_edg))
#print(f'{len(x)} and {len(x[maskCube])}')
x=x[maskCube]
y=y[maskCube]
z=z[maskCube]
rho=rho[maskCube]
h=h[maskCube]
v=v[maskCube]
P=P[maskCube]
B=B[maskCube]
divB=divB[maskCube]
curlB=curlB[maskCube]
mon_est_B=mon_est_B[maskCube]
sol_est_B=sol_est_B[maskCube]
#divBest=divBest[maskCube]

#print(np.max(np.abs(sol_est_B[:,1])),len(sol_est_B[:,1]))
# plot everythinplots(3, 3, sharex=True, figsize=(10, 9))
# plot everything
fig, ax = pl.subplots(2, 3, sharex=True, figsize=(12, 8))
mms=2
aalpha=0.3


n_th=2
widb=0.05
Bx_th=np.exp(-x**2/(2*widb**2))#np.sin(n_th*2*3.1415*x)
divB_th=-2*x/(2*widb**2)*Bx_th#2*3.1415*n_th*np.cos(n_th*2*3.1415*x)


                                                                                                
ax[0][2].plot(x, divB, ".",markersize=mms,alpha=aalpha,label="divB")
ax[0][2].plot(x, curlB[:,2], ".",markersize=mms,alpha=aalpha,label="rotB_z")
ax[0][2].plot(x, divB_th, ".",color='black',markersize=mms,alpha=aalpha,label="divB_th")
ax[0][0].plot(x, Bx_th-np.mean(Bx_th), ".",color='black',markersize=mms,alpha=aalpha,label="B_x_th")
ax[1][2].plot(x,np.abs(mon_est_B[:,0]-np.mean(mon_est_B[:,0]))/np.abs(Bx_th-np.mean(Bx_th)),".",markersize=mms,alpha=aalpha,label="B_mon_est_x/B_x_th")
ax[1][2].plot(x,np.abs(divB)/np.abs(divB_th), ".",markersize=mms,alpha=aalpha,label="divB/divB_th")
    
ax[0][0].plot(x, B[:, 0]-np.mean(B[:,0]), ".", label="d B_x SWIFT",markersize=mms,alpha=aalpha)
ax[0][0].plot(x,mon_est_B[:,0]-np.mean(mon_est_B[:,0]),".",markersize=mms,alpha=aalpha,label="d B_mon_est_x")
ax[0][0].plot(x,sol_est_B[:,0]-np.mean(sol_est_B[:,0]),".",markersize=mms,alpha=aalpha,label="d B_sol_est_x")
ax[0][1].plot(y, B[:, 0]-np.mean(B[:,0]), ".", label="d B_x SWIFT",markersize=mms,alpha=aalpha)
ax[0][1].plot(y,mon_est_B[:,0]-np.mean(mon_est_B[:,0]),".",markersize=mms,alpha=aalpha,label="d B_mon_est_x")
ax[0][1].plot(y,sol_est_B[:,0]-np.mean(sol_est_B[:,0]),".",markersize=mms,alpha=aalpha,label="d B_sol_est_x|")
ax[1][0].plot(x, B[:, 1]-np.mean(B[:,1]), ".", label="B_y-<B_y> SWIFT", markersize=mms,alpha=aalpha)
ax[1][0].plot(x,mon_est_B[:,1]-np.mean(mon_est_B[:,1]),".",markersize=mms,alpha=aalpha,label="d B_mon_est_y")
ax[1][0].plot(x,sol_est_B[:,1]-np.mean(sol_est_B[:,1]),".",markersize=mms,alpha=aalpha,label="d B_sol_est_y")
ax[1][1].plot(y, B[:, 1]-np.mean(B[:,1]), ".", label="B_y-<B_y> SWIFT", markersize=mms,alpha=aalpha)
ax[1][1].plot(y,mon_est_B[:,1]-np.mean(mon_est_B[:,1]),".",markersize=mms,alpha=aalpha,label="d B_mon_est_y")
ax[1][1].plot(y,sol_est_B[:,1]-np.mean(sol_est_B[:,1]),".",markersize=mms,alpha=aalpha,label="d B_sol_est_y")


ax[0][0].set_ylabel("B_x(x)")
ax[0][1].set_ylabel("B_x(y)")
ax[1][0].set_ylabel("B_y(x)")
ax[1][1].set_ylabel("B_y(y)")

mscale_legend=4.0
ax[0][0].legend(loc="best",markerscale=4.0)
ax[0][1].legend(loc="best",markerscale=4.0)
ax[1][0].legend(loc="best",markerscale=4.0)
ax[1][1].legend(loc="best",markerscale=4.0)
ax[0][2].legend(loc="best",markerscale=4.0)
ax[1][2].legend(loc="best",markerscale=4.0)

# Information -------------------------------------
#plt.subplot(236, frameon=False)

text_fontsize = 7

#ax[2][2].text(-0.45, 0.8,
#ax[2][2].plot([-0.45, 0.45], [0.62, 0.62], "k-", lw=1)
#ax[2][2].text(-0.45, 0.5, "$SWIFT$ %s" % git.decode("utf-8"), fontsize=text_fontsize)
#ax[2][2].text(-0.45, 0.4, "$Branch$ %s" % gitBranch.decode("utf-8"), fontsize=text_fontsize)
#ax[2][2].text(-0.45, 0.3, scheme.decode("utf-8"), fontsize=text_fontsize)
#ax[2][2].text(-0.45, 0.2, kernel.decode("utf-8"), fontsize=text_fontsize)
#ax[2][2].text(-0.45, 0.1, "$%.2f$ neighbours" % (neighbours),fontsize=text_fontsize)
#ax[2][2].plot([-0.45, 0.45], [0.0, 0.0], "k-", lw=1)
#ax[2][2].text(-0.45, -0.1, "$\\mu_0:%.2f$ " % (mu0),fontsize=text_fontsize)
#ax[2][2].text(-0.45, -0.2, "$Difussion_\\eta:%.2f$ " % (mhdeta),fontsize=text_fontsize)
#ax[2][2].text(-0.45, -0.3, "$Dedner: Hyp=%.2f // Par=%.2f$" % (dedhyp,dedpar),fontsize=text_fontsize)
#ax[2][2].plot([-0.45, 0.45], [-0.4, -0.4], "k-", lw=1)
#ax[2][2].text(-0.45, -0.5, "$Flavour: $ %s" % mhdflavour.decode("utf-8")[0:30], fontsize=text_fontsize)


# make sure all velocity and magnetic fields plots use the same y axis
# this puts the noise on constant values into context
#for a in [ax[1][0], ax[1][2]]:
#    a.set_ylim(*ax[1][1].get_ylim())
#for a in [ax[2][0], ax[2][2]]:
#    a.set_ylim(*ax[2][1].get_ylim())
ax[0][0].set_ylim(-1.1,1.1)
ax[0][1].set_ylim(-1.1,1.1)
ax[1][0].set_ylim(-1.1,1.1)
ax[1][1].set_ylim(-1.1,1.1)
ax[1][2].set_ylim(0.5,1.1)
# add the time as a title 
ax[0][1].set_title(f"t={t:.2e}")

# save the figure
pl.tight_layout()
pl.savefig(args.output, dpi=300)
