#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Homework 4
# Center of Mass Position and Velocity
# write your name here


# In[2]:


# import modules
import numpy as np
import astropy.units as u
import astropy.table as tbl

from ReadFile import Read


# In[3]:


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot
    
    
    def __init__(self, filename, ptype):
    # Initialize the instance of this Class with the following properties:
    
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        # write your own code to complete this for positions and velocities
        #
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]


    def COMdefine(self,a,b,c,m):
    # Function to compute the center of mass position or velocity generically
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)

        # write your own code to compute the generic COM using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        Acom = np.sum(a * m)/np.sum(m)
        # ycomponent Center of mass
        Bcom = np.sum(b * m)/np.sum(m)
        # zcomponent Center of mass
        Ccom = np.sum(c * m)/np.sum(m)
        
        return Acom, Bcom, Ccom
    
    
    def COM_P(self, delta, VolDec):
    # Function to specifically return the center of mass position and velocity                                         
    # input:                                                                                                           
    #        particle type (1,2,3)                                                                                     
    #        delta (tolerance)                                                                                         
    # returns: One vector, with rows indicating:                                                                                                                                                                            
    #       3D coordinates of the center of mass position (kpc)                                                             

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        # write your own code below
        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # Calculate a new guess for the center of mass by subtracting 
        # the initial guess from the position array
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        RMAX = max(RNEW)/VolDec
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
            # select all particles within the reduced radius (starting from original x,y,z, m)
            #Select values for positions and mass within a set radius
            index2 = np.where(RNEW<RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # Recalculate the COM using the new updated values
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2, y2, z2, m2)
            # Calculate the 
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
            # print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/VolDec
            # check this.                                                                                              
            #print ("maxR", maxR)                                                                                      

            # Calculate a new guess for the center of mass by subtracting 
            # the initial guess from the position array
            xNew = self.x - XCOM2
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                                                       
            COMP = [XCOM, YCOM, ZCOM]

        # set the correct units usint astropy and round all values
        # and then return the COM positon vector
        COMP = np.around(COMP*u.kpc, 2)
        return COMP
    

    def COM_V(self, XCOM,YCOM,ZCOM):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position
        # Calculate a new guess for the center of mass by subtracting 
        # the initial guess from the position array
        xNew = self.x*u.kpc - XCOM
        yNew = self.y*u.kpc - YCOM
        zNew = self.z*u.kpc - ZCOM
        RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)

        
        # Set the new index for particles within the given radius
        indexV = np.where(RNEW<RVMAX)

        # Get the values for the directional velocities within a given radius
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew =  self.m[indexV]
        
        
        # Recalculate the COM Velocity using the new updated values
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew, vynew, vznew, mnew)

        # create a vector to store the COM velocity
        # set the correct units usint astropy
        # round all values
        # Create an array containing the COM velocity values in each direction and then 
        # round the value to 2 decimal places and add units of km/s
        COMV = np.array([VXCOM, VYCOM, VZCOM])
        COMV = np.around(COMV*u.km/u.s,2)
        # return the COM vector                                                                                        
        return COMV


# In[ ]:




