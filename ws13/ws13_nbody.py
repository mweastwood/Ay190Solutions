#!/usr/bin/env python

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as mpl3d

# global constants
ggrav = 6.67e-8
msun  = 1.99e33
seconds_per_year = 24.*3600*365 # roughly
cm_per_pc = 3.1e18
distance_to_sgrAstar = 8e3 * cm_per_pc

# system parameters
id_type = "sun-earth"
if id_type == "sun-earth":
    initial_data_file = "sun_earth.asc"
    distance_unit_to_cm = 1.
    time_unit_to_s = 1.
    mass_unit_to_g = 1.
    Nsteps = 1e4
    t0 = 0
    t1 = 1 * seconds_per_year
    dt = (t1-t0)/Nsteps
elif id_type == "S-stars":
    initial_data_file = "sgrAstar.asc"
    distance_unit_to_cm = np.pi/(180*3600) * distance_to_sgrAstar
    time_unit_to_s = 3600*24*365
    mass_unit_to_g = msun
    Nsteps = 3e4
    t0 = 0
    t1 = 150 * seconds_per_year
    dt = (t1-t0)/Nsteps
else:
    sys.write(sys.stderr, "Unknown id_type %s\n" % id_type)
    sys.exit(1)
final_data_file = "final_positions.asc"

def NbodyRHS(u,mass,time):
    v = np.zeros(u.shape)
    v[:,0:3] = u[:,3:]
    for i in range(0,len(mass)):
        deltaxyz = u[i,0:3] - u[:,0:3]
        r = np.sqrt(np.sum(deltaxyz*deltaxyz,axis=1))
        r[i] = -1 # avoids divivision by zero
        v[i,3] = -ggrav * np.sum(deltaxyz[:,0] * mass/(r**3))
        v[i,4] = -ggrav * np.sum(deltaxyz[:,1] * mass/(r**3))
        v[i,5] = -ggrav * np.sum(deltaxyz[:,2] * mass/(r**3))
    return v

def NbodyRK4(u,mass,time,dt):
    k1 = NbodyRHS(u,mass,time)
    k2 = NbodyRHS(u+0.5*dt*k1,mass,time)
    k3 = NbodyRHS(u+0.5*dt*k2,mass,time)
    k4 = NbodyRHS(u+dt*k3,mass,time)
    return u+dt/6.*(k1+2*k2+2*k3+k4)

def TotalEnergy(u,mass,time):
    (x,y,z,vx,vy,vz) = u.transpose()
    v2 = vx**2+vy**2+vz**2
    Ekin = 0.5*np.sum(mass*v2)
    Egrav = 0.
    for i in range(0,len(mass)):
        deltax = x[i] - x
        deltay = y[i] - y
        deltaz = z[i] - z
        r = np.sqrt(deltax**2 + deltay**2 + deltaz**2)
        r[i] = 1e300 # avoids divivision by zero
        Egrav += - ggrav * mass[i] * np.sum(mass/r)
    return Ekin + Egrav

# main program
plt.ion()

(x,y,z,vx,vy,vz,mass) = np.loadtxt(initial_data_file, unpack = True)


# convert from unitis in initial data file to cgs
x *= distance_unit_to_cm
y *= distance_unit_to_cm
z *= distance_unit_to_cm
vx *= distance_unit_to_cm / time_unit_to_s
vy *= distance_unit_to_cm / time_unit_to_s
vz *= distance_unit_to_cm / time_unit_to_s
mass *= mass_unit_to_g

xmin = np.amin(x)
xmax = np.amax(x)
ymin = np.amin(y)
ymax = np.amax(y)
zmin = np.amin(z)
zmax = np.amax(z)
rmax = 2.5*max(abs(xmin),abs(xmax),abs(ymin),abs(ymax),abs(zmin),abs(zmax))

# use a single state vector to simplify the ODE code
# indices:
# u[:,0] = x
# u[:,1] = y
# u[:,2] = z
# u[:,3] = vx
# u[:,4] = vy
# u[:,5] = vz
u = np.array((x,y,z,vx,vy,vz)).transpose()

for it in range(0, int(Nsteps)):
    time = t0 + it * dt
    u = NbodyRK4(u,mass,time,dt)
    if it % max(1,Nsteps/100) == 0:
      print "it = %d, time = %g years, energy = %g" % \
            (it, time / seconds_per_year,
             TotalEnergy(u,mass,time))
      plt.clf()
      fig = plt.gcf()
      ax = mpl3d.Axes3D(fig)
      ax.scatter(u[:,0],u[:,1],u[:,2])
      ax.set_xlim((-rmax,rmax))
      ax.set_ylim((-rmax,rmax))
      ax.set_zlim((-rmax,rmax))
      plt.draw()

# output result
file_header = "1:x 2:y 3:z 4:vx 5:vy 6:vz 7:mass"
np.savetxt(final_data_file, u, header=file_header)
