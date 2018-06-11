"""
read_migration.py

It reads a csv file that contains the trajectories of several particles
from 3D coordinates and calculates their speed.

Authors:
Jose Guzman, sjm.guzman@gmail.com
Joshua Bagley, sunanjay@imba.oeaw.ac.at
Bajay Sunanjay, joshua.bagley@imba.oeaw.ac.at

Last change: Mon Jun 11 17:58:12 CEST 2018
"""

from __future__ import division

import numpy as np
import pandas as pd

def distance(myarray):
    """
    Calculate the distance between 2 3D coordinates along the 
    axis of the numpy array.
    
    """
    # slice() method is useful for large arrays
    # as in ./local/lib/python2.7/site-packages/numpy/lib/function_base.py
    a = np.asanyarray(myarray)

    slice1 = [slice(None)] # create a slice type object
    slice2 = [slice(None)]
    slice1[-1] = slice(1, None)  # like array[1:] 
    slice2[-1] = slice(None, -1) # like array[:-1]
    slice1 = tuple(slice1)
    slice2 = tuple(slice2)

    # calculate sqrt( dx^2 + dy^2 + dz^2)
    sum_squared = np.sum( np.power(a[slice2]-a[slice1],2), axis=1)
    return np.sqrt( sum_squared )

class Particle(object):
    """
    Uses a pandas object to creates an 'Particle ' object 
    with the speeds and times.

    Usage:
    ------
    >>> import pandas as pd
    >>> mydf = pd.read_csv('.data/myfile.csv', skiprows=2)
    >>> myID = mydf['TrackID'][0] # get first ID
    >>> myparticle = Particle(df = mydf, ID = myID)
    >>> plt.plot(myparticle.time, myparticle.speed)
    """
    def __init__(self, df, ID):
        """
        Reads a custom pandas file

        Arguments:
        ----------
        df -- (panda) is a panda object with locations and TrackIDs
        ID -- (str) the id of the particle to trace
        """
        self.ID = ID
        # select X,Y,Z locations, sampling points and TrackIDs
        sel = df[['Position X','Position Y','Position Z','Time','TrackID']]
        mydf = sel.loc[ sel['TrackID'] == ID ] # row selection

        # Time vector
        self.t = mydf['Time'].values
        self.total_time = np.sum( np.diff(self.t,1) ) * 10 # in minutes
    
        # speed vector
        p = zip(mydf['Position X'], mydf['Position Y'], mydf['Position Z'])

        self.distance = distance(myarray = p) # in um
        self.dt = np.diff(self.t, 1) # time linear interpolation
        self.speed = self.distance/(self.dt*10.) # in um/min

        dt = 1/6. # sampling interval in hours
        self.time = self.t*dt # convert in hours
        self.time = self.time[:-1]# remove last sampling points

        self.samples = len( mydf )

    def __repr__(self):
        """
        Returns readable information of the object by simpy typing it
        """
        info = 'Track ID: %s\n'%self.ID #last sample
        info += 'Traveled time %2.4f min\n'%self.total_time #last sample
        info +='Traveled distance %2.4f um\n'%self.distance.sum() 
        info +='Average speed %2.4f um/min\n'%self.speed.mean() 
        return(info)
        
    def __len__(self):
        """
        Returns the number of samples in the dataset
        """
        return self.samples
        

class CVSReader(object):
    """
    Reads a cvs file containing 3D coordinates to returns particle objects
    """
    def __init__(self, filename=None):
        """
        An object to obtain speed parameters from different particles in 
        in a csv file. It will create a return Particles objects.

        Arguments:
        ----------
        filename -- (str) the path containing the csv file to read

        Usage:
        ------

        >>> from read_migration import reader
        >>> reader('./data/)
        >>> len(reader) # number of particles (e.g. 3)
        >>> a, b, c = reader()
        >>> a.speed # see Particle objects for examples
        """
        self.filename = list() 
        self.particles = list()

        if filename is not None:
            self.filename.append(filename)
            mydf = pd.read_csv(filename, skiprows=2)
            ids = mydf['TrackID'].unique()
            self.particles = [Particle(df = mydf, ID=i) for i in ids]

    def __call__(self, filename=None):
        """
        Return Particles objects form the file

        Arguments:
        ----------
        filename -- (str) the path containing the csv file to read
        """
        if filename is None:
            return self.particles
        else:
            # update attributes
            self.filename = filename
            mydf = pd.read_csv(filename, skiprows=2)
            ids = mydf['TrackID'].unique()
            self.particles = [Particle(df = mydf, ID=i) for i in ids]
            return self.particles

    def __add__(self, filename):
        """
        add two CVSReader objects
        """
        myreader = CVSReader()
        myreader.filename.append(self.filename)

        # particles attribute
        myreader.particles = self.particles
        mydf = pd.read_csv(filename, skiprows=2)
        ids = mydf['TrackID'].unique()
        myreader.particles += [Particle(df = mydf, ID=i) for i in ids]

        return(myreader)
        

    def __len__(self):
        """
        Returns the number of Particle objects 
    
        """
        if self.filename is None:
            return 0
        else:
            return len(self.particles)
            
reader = CVSReader(filename = None)
        

