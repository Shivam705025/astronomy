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
import matplotlib.pyplot as pl
from matplotlib import ticker, cm
matplotlib.use("Agg")
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
    #print(handle["PartType0"].keys())
    gamma = handle["HydroScheme"].attrs["Adiabatic index"][0]
    boxsize = handle["Header"].attrs["BoxSize"][0]
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
    Fmag = handle["PartType0/Fmag"][:]
    Fhyd = handle["PartType0/Fhyd"][:]
    phi = handle["PartType0/DednerScalar"][:]
    git = handle["Code"].attrs["Git Revision"]
    gitBranch = handle["Code"].attrs["Git Branch"]
    scheme = handle["/HydroScheme"].attrs["Scheme"]
    kernel = handle["/HydroScheme"].attrs["Kernel function"]
    neighbours = handle["/HydroScheme"].attrs["Kernel target N_ngb"]
    mu0 = handle["/PhysicalConstants/InternalUnits"].attrs["vacuum_permeability"]
    dedhyp = handle["/HydroScheme"].attrs["Dedner Hyperbolic Constant"]
    dedpar = handle["/HydroScheme"].attrs["Dedner Parabolic Constant"]
    mhdeta = handle["/HydroScheme"].attrs["Diffusion Eta"]
    mhdflavour = handle["/HydroScheme"].attrs["MHD Flavour"]

bb = np.sqrt(B[:,0]*B[:,0]+B[:,1]*B[:,1]+B[:,2]*B[:,2])
vv = np.sqrt(v[:,0]*v[:,0]+v[:,1]*v[:,1]+v[:,2]*v[:,2])
x -= boxsize * 0.5

reg=1e-9
reg_err=1e-2

R0 = np.abs(divB) * h / (bb+reg)
#Fmag = B
#Fmag[:,0] = B[:,0]*divB+(curlB[:,1]*B[:,2]-curlB[:,2]*B[:,1])
#Fmag[:,1] = B[:,1]*divB+(curlB[:,2]*B[:,0]-curlB[:,0]*B[:,2])
#Fmag[:,2] = B[:,2]*divB+(curlB[:,0]*B[:,1]-curlB[:,1]*B[:,0])
ff = np.sqrt(Fmag[:,0]*Fmag[:,0]+Fmag[:,1]*Fmag[:,1]+Fmag[:,2]*Fmag[:,2])
R1 = np.abs((Fmag[:,0]*B[:,0]+Fmag[:,1]*B[:,1]+Fmag[:,2]*B[:,2])/(ff*bb+reg))
cbcb= np.sqrt(curlB[:,0]*curlB[:,0]+curlB[:,1]*curlB[:,1]+curlB[:,2]*curlB[:,2])
R2 = np.abs(divB/(cbcb+reg))

Ftot=Fhyd

fftot=np.sqrt(Ftot[:,0]*Ftot[:,0]+Ftot[:,1]*Ftot[:,1]+Ftot[:,2]*Ftot[:,2])
Rmag=ff/(fftot+ff)

#print(f'the errors are:R0={round(np.mean(R0),2)} R1={round(np.mean(R1),2)} R2={round(np.mean(R2),2)}')
# plot everything
# plot everything
fig, ax = pl.subplots(3, 3, sharex=True, figsize=(18, 14))
mms=2
aalpha=0.3

#change quantity to plot here

def make_color_levels(Q,cmin,cmax,c_res=10,log_sc=True):
    if log_sc:
        levmin=int(np.floor(np.log10(cmin)))
        levmax=int(np.ceil(np.log10(cmax)))
        levels=[]
        levels_short=[]
        for i in range(levmin,levmax):
            levels_short+=[10**i]
            for j in range(int(c_res/10),c_res):
                levels+=[(10/c_res*j)*10**i]
                
    else:
        levels=[cmin+(cmax-cmin)/c_res*k for k in range(c_res)]
        levels_short=[cmin+(cmax-cmin)/c_res*k for k in range(0,c_res,10)]
    return levels, levels_short

def make_density_plot(Q,cmin,cmax,i,j,Q_name,c_res=10,log_sc=True):   
    levels,levels_short=make_color_levels(Q,cmin,cmax,c_res,log_sc)
    if log_sc:
        to_plot=ax[i][j].tricontourf(x,y,Q,levels=np.array(levels),locator=ticker.LogLocator(),cmap = 'viridis')
    else:
        to_plot=ax[i][j].tricontourf(x,y,Q,levels=np.array(levels),cmap = 'viridis')
    fig.colorbar(to_plot,ticks=levels_short)
    ax[i][j].set_ylabel(Q_name)
    return 0

make_density_plot(bb*bb,np.min(bb*bb),1.1*np.max(bb*bb),0,0,"B2",c_res=100,log_sc=False)
#make_density_plot(rho,np.min(rho),np.max(rho),0,1,"rho",c_res=10,log_sc=False)
#make_density_plot(vv,np.min(vv),np.max(vv),0,2,"v",c_res=10,log_sc=False)
make_density_plot(vv,np.min(vv),1.1*np.max(vv),0,2,"|v|",c_res=100,log_sc=False)
#make_density_plot(vv,np.min(vv),1.1*np.max(vv),0,2,"|v|",c_res=100,log_sc=False)
make_density_plot(P,np.min(P),1.1*np.max(P),0,1,"P",c_res=100,log_sc=False)
make_density_plot(np.abs(phi),np.quantile(np.abs(phi),0.001),1.1*np.quantile(np.abs(phi),0.999),1,1,"|phi|",c_res=100,log_sc=False)
#make_density_plot(np.abs(phi),np.quantile(np.abs(phi),0.01),np.quantile(np.abs(phi),0.99),1,1,"|phi|")
make_density_plot(R0+reg_err,1e-2,1,2,0,"R0",c_res=100)
make_density_plot(R1+reg_err,1e-2,1,2,1,"R1",c_res=100)
make_density_plot(R2+reg_err,1e-2,1,2,2,"R2",c_res=100)
#make_density_plot(Rmag+reg_err,1e-2,1,0,2,"Rmag",c_res=100)
make_density_plot(Rmag+reg_err,1e-2,1,1,0,"R_mag",c_res=100,log_sc=False)
#make_density_plot(R0,np.quantile(R0,0.001),1.1*np.quantile(R0,0.999),2,0,"R0")
#make_density_plot(R1,np.quantile(R1,0.001),1.1*np.quantile(R1,0.999),2,1,"R1")
#make_density_plot(R2,np.quantile(R2,0.001),1.1*np.quantile(R2,0.999),2,2,"R2")




text_fontsize=7

ax[1][2].text(-0.45, 0.8,
    "OrszagTangVortex with $\\gamma=%.3f$ in 3D\nat $t=%.2f$" % (gamma, t),
    fontsize=text_fontsize,
)
ax[1][2].text(-0.45, 0.75, "$SWIFT$ %s" % git.decode("utf-8"), fontsize=text_fontsize)
ax[1][2].text(-0.45, 0.7, "$Branch$ %s" % gitBranch.decode("utf-8"), fontsize=text_fontsize)
ax[1][2].text(-0.45, 0.65, scheme.decode("utf-8"), fontsize=text_fontsize)
ax[1][2].text(-0.45, 0.6, kernel.decode("utf-8"), fontsize=text_fontsize)
ax[1][2].text(-0.45, 0.55, "$%.2f$ neighbours" % (neighbours),fontsize=text_fontsize)
ax[1][2].text(-0.45, 0.5, "$%.0f$ particles" % (len(x)),fontsize=text_fontsize)
ax[1][2].text(-0.45, 0.45, "$\\mu_0:%.2f$ " % (mu0),fontsize=text_fontsize)
ax[1][2].text(-0.45, 0.4, "$Difussion_\\eta:%.2f$ " % (mhdeta),fontsize=text_fontsize)
ax[1][2].text(-0.45, 0.35, "$Dedner: Hyp=%.2f // Par=%.2f$" % (dedhyp,dedpar),fontsize=text_fontsize)
ax[1][2].text(-0.45, 0.3, "$Flavour: $ %s" % mhdflavour.decode("utf-8")[0:30], fontsize=text_fontsize)
ax[1][2].text(-0.45, 0.2, "R0 = h divB/B, R1 = Cos(F_mag^B), R2 = divB / |rotB|, Rmag = |F_mag|/(|F_hydro|+|F_mag|)", fontsize=text_fontsize)
ax[1][2].text(-0.45, 0.15, "<R0> = $%.2f$, <R1> = $%.2f$, <R2> = $%.2f$, <Rmag> = $%.2f$" % (np.mean(R0),np.mean(R1),np.mean(R2),np.mean(Rmag)), fontsize=text_fontsize)
ax[1][2].text(-0.45, 0.1, "max(R0) = $%.2f$, max(R1) = $%.2f$, max(R2) = $%.2f$, max(Rmag) = $%.2f$" % (np.max(R0),np.max(R1),np.max(R2),np.max(Rmag)), fontsize=text_fontsize)


# add the time as a title
ax[0,1].set_title(f"t={t:.2e}")

# save the figure
pl.tight_layout()
pl.savefig(args.output, dpi=300)
