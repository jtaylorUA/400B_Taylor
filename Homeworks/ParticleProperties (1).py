#!/usr/bin/env python
# coding: utf-8

# In[100]:


import numpy as np
import astropy.units as u
#from ReadFile import Read
def ParticleInfo(filename,particle_type,particle_number):
    file = open(filename, 'r') #Opens the file in read mode
    line1 = file.readline() #Pics the first line in the file and reads it
    label, value = line1.split() #Splits the line at a blank space and stores it as the label and the value
    time = float(value)*u.Myr #Reads the value in units of Myr
    line2 = file.readline() #Reads the second line in the file and splits the line in the same way as the first line
    label, value = line2.split()
    number = float(value)
    #Store the remaining data from the file so that it can be accessed using the type of data and which line 
    data = np.genfromtxt(filename,dtype=None,names = True,skip_header=3) 
    file.close()
    def Get_value(particle_number):
        #Gets the mass of the 100th particle (solMass)
        mass = data['m'][particle_number]
        #Pull the data for positional data: (x,y,z)
        #The values are squared for ease of use in calculation of 3D distance (kpc)
        x = data['x'][particle_number]**2
        y = data['y'][particle_number]**2
        z = data['z'][particle_number]**2
        #Pull the values for the directional velocities (vx,vy,vz) (km/s)
        #The values are squared for ease of use in calculation of 3D velocity
        vx = data['vx'][particle_number]**2
        vy = data['vy'][particle_number]**2
        vz = data['vz'][particle_number]**2
        #Calculate the 3D distance and velocity
        dist = np.sqrt((x+y+z))*u.kpc
        #Converts distance to units of lightyears using .to function in astropy
        dist = dist.to(u.lyr)
        vel = np.sqrt((vx+vy+vz))
        print('3D Distance: ',dist)
        print('3D Velocity: ',vel,'km/s')
        print('Mass: ',mass,'solMass')
    Get_value(99)
ParticleInfo('Downloads/MW_000.txt', 'disk', 101)  


# In[ ]:





# In[ ]:




