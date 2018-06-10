"""
nerst.py

Jose Guzman, sjm.guzman@gmail.com
Created: Fri Apr 29 08:04:57 CEST 2016

Solves the Nerst equation for different chloride conditions.

"""
from math import log
from terminaltables import AsciiTable 
import numpy as np

def nerst(t, oution, inion, z):
    """ 
    Solves the following equation:

    .. math:: 
        E = \frac{R T}{z F} \ln\frac{out}{in}\\

    Arguments
    ---------

    oution : float
             Extracellular ionic concentration
    intion : float 
             Intracellular ionic concentration
    z      : int
             Valence of the ion
    t      : temp
             Temperature in celsius

    Returns
    -------
    voltage : float 
            the potential (in mV) at which the net flux of current is zero.
        

    Examples
    --------
    >>> nerst( t = 30, oution = 126.5, inion = 10, z=-1)
    >>> -66.294937

    """
    K = 273.16 + t # transform in kelvin
    volt = ((8.31451*K)/(z*96485.3))*log(oution/inion) 
    return(volt*1000) # in mV

if __name__ == '__main__':
    data = [['Reference', 'E_rev (Cl-) mV']]
    # data from Pavlidis % Madison 1999
    ECl = nerst(t = 30, oution=126.5, inion=10, z=-1)
    data.append(['Pavlidis & Madison, 1999', ECl])

    # data from Sasaki et al
    ECl = nerst(t = 30, oution=133.3, inion=4, z=-1)
    data.append(['Sasaki et al., 2012', ECl])

    # data from Mitra et al, 2011
    ECl = nerst(t = 20, oution=129, inion=9, z=-1)
    data.append(['Mitra et al., 2011', ECl])

    # data from Kraushaar and Jonas, 2000 
    ECl = nerst(t = 20, oution=126.5, inion=149, z=-1)
    data.append(['Krausharr and Jonas, 2000',ECl])

    # data from Espinoza 
    ECl = nerst(t = 20, oution=133.5, inion=44, z=-1)
    data.append(['Espinoza et al., **', ECl])

    table = AsciiTable(data)
    print table.table

    import matplotlib.pyplot as plt
    x = np.arange(0.1, 50.0, 0.01)
    k = lambda x:58*np.log10(x/100.) # K-nerst equation
    y = k(x)
    plt.semilogx(x,y, color='royalblue')
    plt.vlines(2.5, -120, k(2.5), linestyle=':', color='brown')
    plt.hlines(k(2.5), 0.01, 2.5, linestyle=':', color='brown')
    plt.ylim(ymin=-120)
    plt.xlim(xmin=0.1)

    plt.ylabel('Resting membrane potential (mV)')
    plt.xlabel('Log [K$^+$]')
    plt.show()



