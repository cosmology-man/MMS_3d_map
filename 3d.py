# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 15:34:28 2020

@author: Asha
"""
import math as m
import numpy as np


n = np.array(range(0,32761))


#point on the surface
r = 6371
h = np.array(range(0, 181))
th = np.array(range(-90, 91))
n_thh = range(0, 181)

#surface vector
spherical_array = []
for i in h:
    phi = i
    for i in th:
        vector = [r, phi, i]
        spherical_array.append(vector)

#cartesian
cartesian_array = []
for i in spherical_array:
    x = r*m.sin(m.radians(i[1]))*m.cos(m.radians(i[2]))
    y = r*m.sin(m.radians(i[1]))*m.sin(m.radians(i[2]))
    z = r*m.cos(m.radians(i[1]))
    vector = [x, y, z]
    cartesian_array.append(vector)

#center to sun vector
l = [1496004241, 0, 0]

#surface to sun vector
sts_array = []
for i in cartesian_array:
    x = l[0] - i[0]
    y = l[1] - i[1]
    z = l[2] - i[2]
    vector = [x, y, z]
    sts_array.append(vector)

#sts dotted with surface
adotb_array = []
for i in n:
    x = cartesian_array[i][0] * sts_array[i][0]
    y = cartesian_array[i][1] * sts_array[i][1]
    z = cartesian_array[i][2] * sts_array[i][2]
    adotb = x + y + z
    adotb_array.append(adotb)

#sts length
sts_len = []
for i in sts_array:
    length = m.sqrt(i[1]**2 + i[0]**2 + i[2]**2)
    sts_len.append(length)

#surface length
surf_len = []
for i in cartesian_array:
    length = m.sqrt(i[1]**2 + i[0]**2 + i[2]**2)
    surf_len.append(length)

#theta
theta_array = []
for i in n:
    lal_lbl = surf_len[i]*sts_len[i]
    theta = adotb_array[i]/lal_lbl
    theta = m.acos(theta)
    theta_array.append(theta) 

theta = []
for i in theta_array:
    convert = m.degrees(i)
    theta.append(convert)

#Bz
Bz = 0.012695087
Bz_sig = 3.2014756

#Dp
Dp = 1.7389104
Dp_sig = 1.4139446

#variations
Bzplus = Bz+Bz_sig
Dpplus = Dp+Dp_sig
Bzminus = Bz-Bz_sig
Dpminus = Dp-Dp_sig

#variation arrays
dplu = [Bzplus, Dpplus]
fplu = [Bzplus, Dpminus]
dmin = [Bzminus, Dpminus]
fmin = [Bzminus, Dpplus]
over_var = [dplu, fplu, dmin, fmin]

#r0
r0pp = 10.22+1.29* m.tanh(0.184*(Bzplus+8.14))*(Dpplus**(-1/6.6))
r0p = 10.22+1.29* m.tanh(0.184*(Bzplus+8.14))*(Dpminus**(-1/6.6))
r0m = 10.22+1.29* m.tanh(0.184*(Bzminus+8.14))*(Dpplus**(-1/6.6))
r0mm = 10.22+1.29* m.tanh(0.184*(Bzminus+8.14))*(Dpminus**(-1/6.6))

#alpha
app = (0.58 - 0.007*Bzplus)*(1 + 0.024*m.log(Dpplus))
apm = (0.58 - 0.007*Bzplus)*(1 + 0.024*m.log(Dpminus))
amm = (0.58 - 0.007*Bzminus)*(1 + 0.024*m.log(Dpminus))  
amp = (0.58 - 0.007*Bzminus)*(1 + 0.024*m.log(Dpplus))

#variable set
set1 = [r0pp, app]
set2 = [r0p, apm]
set3 = [r0m, amp]
set4 = [r0mm, amm]

#magnetopause distance
radius_array = []
for i in theta_array:
    r1 = set1[0]*(2/(1+m.cos(i)))**set1[1]
    r2 = set2[0]*(2/(1+m.cos(i)))**set2[1]
    r3 = set3[0]*(2/(1+m.cos(i)))**set3[1]
    r4 = set4[0]*(2/(1+m.cos(i)))**set4[1]
    rt = [r1, r2, r3, r4]
    radius_array.append(rt)

radius = []
for i in radius_array:
    rmin = min(i)
    rmax = max(i)
    rf = [rmin, rmax]
    radius.append(rf)
print(min(radius))
print(max(radius))
    
#unitize directional vector
unit_vector = []
for i in n:
    u1 = cartesian_array[i][0]/surf_len[i]
    u2 = cartesian_array[i][1]/surf_len[i]
    u3 = cartesian_array[i][2]/surf_len[i]
    u = [u1, u2, u3]
    unit_vector.append(u)

graph_coordinates = []
for i in n:
    g1 = unit_vector[i][0]*radius[i][0]
    g2 = unit_vector[i][1]*radius[i][0]
    g3 = unit_vector[i][2]*radius[i][0]
    g = [g1, g2, g3]
    graph_coordinates.append(g)

#graph
#%matplotlib notebook
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection = '3d')

X = []
Y = []
Z = []
for i in graph_coordinates:
    X.append(i[0])
    Y.append(i[1])
    Z.append(i[2])
Axes3D.plot_trisurf(ax, X, Y, Z)


      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      