#!/usr/bin/env python
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
import numpy as np
import h5py
import matplotlib.pyplot as plt
import pickle 
# from scipy.integrate import odeint

t = np.linspace(0, 40, 100000)
y0 = [0, 10]
a = 30.0
G = 4.300927e-06
M = 1e15
GM = G * M

lengthrun = 201
numbpar = 3

radius = np.zeros((numbpar, lengthrun))
xx = np.zeros((numbpar, lengthrun))
yy = np.zeros((numbpar, lengthrun))
zz = np.zeros((numbpar, lengthrun))
time = np.zeros(lengthrun)
for i in range(0, lengthrun):
    Data = h5py.File("output/output_%04d.hdf5" % i, "r")
    header = Data["Header"]
    time[i] = header.attrs["Time"]
    particles = Data["PartType1"]
    positions = particles["Coordinates"]
    xx[:, i] = positions[:, 0] - 500.0
    yy[:, i] = positions[:, 1] - 500.0
    zz[:, i] = positions[:, 2] - 500.0

col = ["b", "r", "c", "y", "k"]

#Plots the orbits
fig,ax = plt.subplots(nrows=1, ncols=3, num=1, figsize=(12, 4.1))
fig.suptitle("Orbits", fontsize=15)
ax[0].clear() ; ax[1].clear()

for i in range(0, numbpar):
    ax[0].plot(xx[i, :], yy[i, :], col[i])

ax[0].set_aspect('equal', 'box')
ax[0].set_xlim([-35,35])
ax[0].set_ylim([-35,35])
ax[0].set_ylabel("y (kpc)")
ax[0].set_xlabel("x (kpc)")


for i in range(0, numbpar):
    ax[1].plot(xx[i, :], zz[i, :], col[i])

ax[1].set_aspect('equal', 'box')
ax[1].set_xlim([-35,35])
ax[1].set_ylim([-35,35])
ax[1].set_ylabel("z (kpc)")
ax[1].set_xlabel("x (kpc)")
ax[1].legend(["Particule 1, $R = 5$ kpc", "Particule 2, $R = 5$ kpc", "Particule 3, $R = 30$ kpc"])

for i in range(0, numbpar):
    ax[2].plot(yy[i, :], zz[i, :], col[i])

ax[2].set_aspect('equal', 'box')
ax[2].set_xlim([-35,35])
ax[2].set_ylim([-35,35])
ax[2].set_ylabel("z (kpc)")
ax[2].set_xlabel("y (kpc)")
plt.savefig("circular_orbits.png")
plt.close()

#%%
#Plot of the deviation from circular orbit
R_1_th = 5.0 #kpc
R_2_th = 5.0 #kpc
R_3_th = 30.0 #kpc

fig2, ax2 = plt.subplots(nrows=1, ncols=2, num=2, figsize=(12, 4.5))
fig2.suptitle("Deviation from circular orbit", fontsize=15)
ax2[0].clear() ; ax2[1].clear()

#Gather the x,y and z components of each particule into one array
pos_1 = np.array([xx[0, :], yy[0, :], zz[0, :]])
pos_2 = np.array([xx[1, :], yy[1, :], zz[1, :]])
pos_3 = np.array([xx[2, :], yy[2, :], zz[2, :]])

#Compute the radii
r_1 = np.linalg.norm(pos_1, axis=0) ; error_1 = np.abs(r_1 - R_1_th)/R_1_th*100
r_2 = np.linalg.norm(pos_2, axis=0) ; error_2 = np.abs(r_2 - R_2_th)/R_2_th*100
r_3 = np.linalg.norm(pos_3, axis=0) ; error_3 = np.abs(r_3 - R_3_th)/R_3_th*100

ax2[0].plot(time, error_1, col[0])
ax2[1].plot(time, error_2, col[1])
ax2[1].plot(time, error_3, col[2])
ax2[0].set_ylabel("Deviation (\%)") ; ax2[0].set_xlabel("Time (Gyr)")
ax2[1].set_ylabel("Deviation (\%)") ; ax2[1].set_xlabel("Time (Gyr)")
ax2[0].legend(["Particule 1, $R = 5$ kpc"])
ax2[1].legend(["Particule 2, $R = 5$ kpc", "Particule 3, $R = 30$ kpc"])

plt.savefig("deviation.png")
plt.close()

#%%Make a comparison with the obtained data and ours to check nothing is broken
#Then delete the code to save the reference data

filename = "original_radii.pkl" #Original data
# with open(filename, "wb") as file:
#     pickle.dump(np.array([r_1, r_2, r_3]), file)
#     # pickle.dump(r_2, file)
#     # pickle.dump(r_3, file)

with open(filename, "rb") as file:
    radii_original = pickle.load(file)
r_1_exp = radii_original[0, :] #exp stands for "experimental"
r_2_exp = radii_original[1, :]
r_3_exp = radii_original[2, :]

fig3, ax3 = plt.subplots(nrows=1, ncols=3, num=3, figsize=(12, 4.3))
fig3.suptitle("Deviation from the original data", fontsize=15)
ax3[0].clear() ; ax3[1].clear() ; ax3[2].clear() 

error_1 = np.abs(r_1 - r_1_exp)/r_1_exp*100 
error_2 = np.abs(r_2 - r_2_exp)/r_2_exp*100 
error_3 = np.abs(r_3 - r_3_exp)/r_3_exp*100 

ax3[0].plot(time, error_1, col[0])
ax3[1].plot(time, error_2, col[1])
ax3[2].plot(time, error_3, col[2])
ax3[0].set_ylabel("Deviation (\%)") ; ax3[0].set_xlabel("Time (Gyr)")
ax3[1].set_ylabel("Deviation (\%)") ; ax3[1].set_xlabel("Time (Gyr)")
ax3[2].set_ylabel("Deviation (\%)") ; ax3[2].set_xlabel("Time (Gyr)")
ax3[0].legend(["Particule 1, $R = 5$ kpc"])
ax3[1].legend(["Particule 2, $R = 5$ kpc"])
ax3[2].legend(["Particule 3, $R = 30$ kpc"])
plt.savefig("deviation_from_original_data.png")
plt.close()
