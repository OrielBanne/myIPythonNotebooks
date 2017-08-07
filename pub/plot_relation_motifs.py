"""
plot_relation_motifs.py

Jose Guzman, sjm.guzman@gmail.com
Date: Sat Jul 22 15:36:05 CEST 2017
Last change: Tue Jul 25 10:19:07 CEST 2017

Plots the relation between synaptic weights of all possible 2-edge
connectivity motifs found between CA3 neurons. For example, for a 
bidirectional connection of the type

    [syn1]
A <--------> B
    [syn2]

A point with coordinates (syn1, syn2) will be plotted. Thus, for the
plot will contain one coordinate for every pair of connections. 

It is based on the data Guzman et al., 2016 and includes:

6 bidirectional motifs  (12 coordinates)
20 divergent motifs (40 coordinates)
10 convergent motifs (20 coordinates)
24 linear motifs ( 18 coordinates)

"""

import numpy as np

from scipy.stats import linregress
from scipy.stats import t as T

import matplotlib.pyplot as plt

from terminaltables import AsciiTable

#--------------------------------------------------------------------------
# Load synaptic data in a list of tuples
#-------------------------------------------------------------------------- 

# 6 bidirectional motifs: 12 synapses
bidirectional = [
    (15.7548, 28.2989), # (51, 53)
    (11.153, 34.7268),  # (82, 83)
    (6.77177, 32.841),  # (91, 95)
    (31.1624, 34.616),  # (92, 96)
    (7.56688, 6.10144), # (128, 131)
    (6.32708, 7.82399)  # (129, 133)
]

# 10 convergent motifs: 20 synapses 
convergent = [ 
    (11.7306, 15.3488), # (67, 68)
    (9.89007, 13.7035), # (70, 71)
    (11.4741, 11.153),  # (81, 82) 
    (6.77177, 3.58013), # (91, 97)
    (32.841, 34.616),   # (95, 96)
    (4.10465, 6.10144), # (126, 130)
    (4.10465, 6.8877),  # (126, 131)
    (6.10144, 17.8282), # (130, 132)
    (7.56688,  6.01477),# (128, 135)
    (6.8877, 7.82399)   # (131, 133)
]

# 20 divergent motifs: 40 synapses
divergent = [
    (13.2576, 15.3028), # (25, 26)
    (11.0663, 28.2989), # (52, 53)
    (8.32979,15.6352),  # (74, 75)
    (6.77177, 31.5107), # (91, 92)
    (6.77177, 31.1624), # (91, 93)
    (6.77177, 8.45281), # (91, 94)
    (31.5107, 31.1624), # (92, 93)
    (31.5107, 8.45281), # (92, 94)
    (31.1624, 8.45281), # (93, 94)
    (20.0786, 4.55908), # (120, 121)
    (4.10465, 29.9351), # (126, 127)
    (4.10465, 7.56688), # (126, 128)
    (4.10465, 6.32708), # (126, 129) 
    (29.9351, 7.56688), # (127, 128)
    (29.9351, 6.32708), # (127, 129)
    (7.56688, 6.32708), # (128, 129)
    (6.8877, 6.10144),  # (131, 130)
    (7.82399, 17.8282), # (133, 132)
    (7.82399, 3.63378), # (133, 134)
    (17.8282, 3.63378)  # (132, 134)
]

#  24 linear motifs: 48 synapses
linear = [
    (5.95899, 16.9644), # (4,5)
    (14.7195, 16.8318), # (20, 21) 
    (15.7548, 11.0663), # (51, 52)
    (11.4741, 34.7268), # (81, 83)
    (6.77177, 31.5107), # (91, 92)
    (6.77177, 31.1624), # (91, 93)
    (6.77177, 8.45281), # (91, 94)
    (34.616, 6.77177),  # (96, 91)
    (34.616, 31.1624),  # (96, 93)
    (34.616, 8.45281),  # (96, 94) 
    (3.58013, 6.77177), # (97, 91)
    (11.396, 12.2366),  # (102, 103)
    (7.56688, 6.10144), # (128, 130)
    (6.32708, 17.8282), # (129, 132)
    (6.32708, 3.63378), # (129, 134)
    (6.8877, 4.10465),  # (131, 126)
    (6.8877, 29.9351),  # (131, 127)
    (6.8877, 6.32708),  # (131, 129)
    (7.82399, 4.10465), # (133, 126)
    (7.82399, 29.9351), # (133, 127)
    (7.82399, 7.56688), # (133, 128)
    (6.01477, 6.10144), # (135, 130)
    (6.01477, 6.8877),  # (135, 131)
    (3.63378, 6.01477)  # (134, 135)
    
]

#-------------------------------------------------------------------------
# Plot data and linear regression together, 
# and output table
#-------------------------------------------------------------------------
def plot_linear_fit(xdata, ydata, color = None, title = None, ax = None):
    """
    Plots the linear fit togheter with the two-side 95% confident intervals 
    
    Parameters
    ----------
    xdata: 1D Numpy array
    ydata: 1D Numpy array 
    color: color object for the linear fit [1], eg = 'red'
    title: (string) containing the title of the output table
    
    Returns
    -------
    A table with the parameters for linear fit of xdata to ydata.

    For 95% confident intervals see:
    https://tomholderness.wordpress.com/2013/01/10/confidence_intervals/

    [1] see https://matplotlib.org/examples/color/named_colors.html
    
    """
    if ax is None:
        ax = plt.gca() # if not give, get current axis
        
    m, a, rval, pval, stderr = linregress(xdata, ydata)

    # the linear function
    f = lambda(x): a + m*x # linear function
    xfit = np.linspace(np.min(xdata), np.max(xdata),100)
    yfit = f(xfit)

    y_err = ydata - f(xdata) # residuals
    SSE = np.power(y_err, 2).sum() # sum of squared errors

    # Calculate confident intervals
    mean_x = np.mean(xdata)	   # mean of data
    n = xdata.size             # number of samples
    # for a 2 tailed 95% confidence interval enter 0.975
    tstat = T.ppf(0.975, n-1)  # appropriate t value

    confs = tstat * np.sqrt( (SSE/(n-2)) * (1.0/n + 
        (np.power((xfit-mean_x),2)/ ((np.sum(np.power(xdata,2)))-n*(np.power(mean_x,2))))))
    
    lower_conf = yfit - abs(confs)
    upper_conf = yfit + abs(confs)

    ax.plot(xfit, yfit, lw =2 , color = color)
    ax.plot(xfit, lower_conf, '--', lw = 1, color = color)
    ax.plot(xfit, upper_conf, '--', lw = 1, color = color)
    ax.fill_between(xfit,upper_conf,lower_conf,  color = color, alpha=.1)

    ax.text(28, 40,'P = %2.4f'%pval, color = color)

    # stdout statistic
    infostats = [
        (title, 'value'),
        ('Slope', m), 
        ('Intercept', a), 
        ('Correlation coef', rval), 
        ('P-value (slope=0)', pval),
        ('Standard error', stderr)
    ]

    print AsciiTable(infostats).table

def plot_linearfit(data, color = None, title = None, ax = None ):
    """
    Fit linearly the data from the pairs of synapses in all possible
    two-edges motifs (e.g. bidirectional, convergent, divergent and
    linear motifs).

    data    -- a tuple containing pairs of synapses in motifs
    color   -- a matplotlib color [1] 
    ax      -- axis object

    Returns:
    an axis object with the data, a linear fit with 95% confident
    interval, and a table with statistics from the linear regression. 

    [1] see https://matplotlib.org/examples/color/named_colors.html
    """
    if ax is None:
        ax = plt.gca() # if not give, get current axis
    
    xdata, ydata = np.array ( zip(*data) )
    plot_linear_fit(xdata, ydata, color, title, ax)

    ax.plot(xdata, ydata, 'o', markersize = 5, color = color)
    
    ax.set_title(title, color = color)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis='both', direction='out')
    ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
    ax.get_yaxis().tick_left()

    ax.set_xlim(-1, 45)
    ax.set_ylim(-1,45)
    ax.set_xlabel('Syn1 (pA)')
    ax.set_ylabel('Syn2 (pA)')

    return(ax)

# plot linear regression of four motif types
fig, ax2D = plt.subplots(2,2)
plot_linearfit(data = bidirectional, 
    color = 'brown', title='Bidirectional motifs', ax = ax2D[0,0])

plot_linearfit(data = convergent, 
    color = 'darkgreen', title='Convergent motifs', ax = ax2D[0,1])

plot_linearfit(data = divergent, 
    color = 'royalblue', title='Divergent motifs', ax = ax2D[1,0])

plot_linearfit(data = linear, 
    color = 'purple', title='Linear motifs', ax = ax2D[1,1])


# Fine-tune figure 
# hide x ticks and label for top plots 
plt.setp([axis.get_xticklabels() for axis in ax2D[0, :]], visible=False)

# hide  y ticks for right plots
plt.setp([axis.get_yticklabels() for axis in ax2D[:, 1]], visible=False)

plt.show()
