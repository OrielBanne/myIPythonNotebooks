"""
migration.py

It reads a custom csv file that contains the trajectories of several 
particles in 3D coordinates. It calculates the speed and their
trajectories.

Authors:
Jose Guzman, sjm.guzman@gmail.com
Joshua Bagley, sunanjay@imba.oeaw.ac.at
Bajay Sunanjay, joshua.bagley@imba.oeaw.ac.at

Last change: Tue Jun 12 05:07:19 CEST 2018
"""

from __future__ import division

import numpy as np
import pandas as pd

def distance(myarray):
    """
    Calculate the distance between two 3D coordinates along the 
    first axis of the numpy array.

    Arguments:
    ----------
    myarray (NumPY array)
    
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
    Uses a pandas dataframe to creates a 'Particle' object 
    that contains the speeds and times.

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
        Reads a custom pandas dataframe containing  
        'Position X','Position Y','Position Z','Time' and 'TrackID'.

        Arguments:
        ----------
        df -- (panda) is a panda dataframe with locations and TrackIDs
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
        

class ParticleReader(object):
    """
    Reads a cvs file containing 3D coordinates to manage Particle objects
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
        # without filename attributes are empty lists
        self.filename = list() 
        self.particle = list()

        if filename is not None:
            self.filename.append(filename)

            mydf = pd.read_csv(filename, skiprows=2)
            ids = mydf['TrackID'].unique()
            self.particle = [Particle(df = mydf, ID=i) for i in ids]

    def __call__(self, filename=None):
        """
        Return Particles objects form the file

        Arguments:
        ----------
        filename -- (str) the path containing the csv file to read

        Returns:
        -------
        A ParticleReader Object
        """

        return ParticleReader( filename )

    def __add__(self, ParticleReaderObj):
        """
        add two ParticleReader objects

        Arguments:
        ----------
        ParticleReaderObj -- A ParticleReaderObj
        """

        myObj = ParticleReaderObj 
        myreader = ParticleReader(None) 

        if not myObj.filename: # reader empty 
            print('empty object')
            return self

        # add filenames
        myreader.filename = self.filename + myObj.filename

        # add Particle objects
        myreader.particle = self.particle +  myObj.particle

        return( myreader )
        

    def __len__(self):
        """
        Returns the number of Particle objects 
    
        """
        return len(self.particle)
            
reader = ParticleReader(filename = None)
        

