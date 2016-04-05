# Copyright 2016 Jochen Kursawe. See the LICENSE file at the top-level directory 
# of this distribution and at https://github.com/kursawe/MCSTracker/blob/master/LICENSE.

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from os import path
from os.path import dirname
from scipy.optimize import curve_fit
from nose.plugins.attrib import attr

"""Obtain averages for the times it takes to run mcs finding"""

def power_law(x, a, b):
    return a * np.power(x,b) 

mpl.rcParams['mathtext.default'] = 'regular'
# make a dictionary with mesh sizes as keys and computation times as Values

sizes_and_times = np.loadtxt('parallel_timings_array.csv')

# sizes_and_times = sizes_and_times[sizes_and_times[:,0].argsort()]

axis_lengths = sizes_and_times[:,0]
calculation_times = sizes_and_times[:,1:]

sizes = [length*length for length in axis_lengths]

popt, pcov = curve_fit(power_law, sizes, np.mean(calculation_times,axis =1))

print 'in a * x^b a is'
print popt[0]
print 'and b is'
print popt[1]

figuresize = (4,2.75)
font = {'size'   : 10}
plt.rc('font', **font)

print len(calculation_times[0,:])
print np.shape(calculation_times)
x_values = np.random.normal(0, 4, size = len(calculation_times[0,:]))
print len(x_values)
this_alpha = 0.7

fine_x_values = np.linspace(0, max(sizes), 100)

plot_figure = plt.figure(figsize = figuresize)
plt.xlabel('Number of cells in tissue n')
plt.ylabel('Calculation time / seconds')
data_markers = plt.errorbar(sizes, np.mean(calculation_times,axis = 1), 
             yerr = np.std(calculation_times,axis = 1), fmt = '.')
for index, size in enumerate(sizes):
    plt.plot( size + x_values, calculation_times[index], color='orange',marker='.', alpha = this_alpha,lw=0, zorder=1 )
fit_marker, = plt.plot(fine_x_values, power_law(fine_x_values, popt[0], popt[1]) )
plt.xlim(0,1700)
plt.locator_params(axis = 'x', nbins=6)   
plt.locator_params(axis = 'y', nbins=4)   
plt.legend( [data_markers, fit_marker], ('Measured', 'Fit $a\cdot n^b$'), loc = 'upper left', fontsize = 10)
plot_figure.tight_layout()
file_name = "calculation_times.pdf"
plot_figure.savefig(file_name)

# log_figure = plt.figure(figsize = figuresize)
# plt.xlabel('Mesh size')
# plt.ylabel('Calculation time')
# plt.semilogy(sizes_and_times[:,0], sizes_and_times[:,1])
# plt.semilogy(sizes_and_times[:,0], power_law(sizes_and_times[:,0], popt[0], popt[1]) )
# print sizes_and_times
# log_figure.tight_layout()
# file_name = "log_times.pdf"
# log_figure.savefig(path.join(dirname(__file__),'output',file_name))
# 
# log_log_figure = plt.figure(figsize = figuresize)
# plt.xlabel('Mesh size')
# plt.ylabel('Calculation time')
# plt.loglog(sizes_and_times[:,0], sizes_and_times[:,1])
# plt.loglog(sizes_and_times[:,0], power_law(sizes_and_times[:,0], popt[0], popt[1]) )
# print sizes_and_times
# log_log_figure.tight_layout()
# file_name = "log_log_times.pdf"
# log_log_figure.savefig(path.join(dirname(__file__),'output',file_name))
