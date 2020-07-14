# Same as pdot.py but we work w/ the orbital separation.
# Differentiating Kepler 3rd leads to :
# (adot/a) / (-Mdot/M) = (2/3) * (Pdot/P) / (-Mdot/M) - a * q / (1+q)
# where M is the mass of the donor star
# N.B. : we keep writing Pdot in the script below to minimize the changes

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.cm as cm
from matplotlib.colors import rgb2hex
# from bokeh import mpl
# from bokeh.plotting import output_file, show, ColumnDataSource, figure, vplot
# from bokeh.models import HoverTool
import pandas as pd
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FormatStrFormatter
import matplotlib.ticker as ticker
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
import holoviews as hv
# import colorcet as cc
from matplotlib.colors import LinearSegmentedColormap

from matplotlib.patches import Rectangle

import os
import re

import subprocess

from parfile import *

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plots the CEDE median factor (in radial log bins) averaged over the 2nd half of the simulation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

hv.notebook_extension("matplotlib")

def fmt(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	return r'{} $\cdot $10$^{{{}}}$'.format(a, b)

def fmtCbar(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	if (b!=0 and b!=1 and b!=-1):
		return r'{}E{}'.format(a, b)
	elif (b==0):
		return r'{}'.format(a)
	elif (b==1):
		return '{:2.0f}'.format(x)
	elif (b==-1):
		return '{:3.2f}'.format(x)

font = {'family' : 'normal',
'size'   : fontsize}
#'weight' : 'bold',

plt.rc('font', **font)

def add_arrow(line, position=None, direction='right', size=15, color=None):
    """
    add an arrow to a line.

    line:       Line2D object
    position:   x-position of the arrow. If None, mean of xdata is taken
    direction:  'left' or 'right'
    size:       size of the arrow in fontsize points
    color:      if None, line color is taken.
    """
    if color is None:
        color = line.get_color()

    xdata = line.get_xdata()
    ydata = line.get_ydata()

    if position is None:
        position = xdata.mean()
    # find closest index
    start_ind = np.argmin(np.absolute(xdata - position))
    if direction == 'right':
        end_ind = start_ind + 1
    else:
        end_ind = start_ind - 1

    line.axes.annotate('',
        xytext=(xdata[start_ind], ydata[start_ind]),
        xy=(xdata[end_ind], ydata[end_ind]),
        arrowprops=dict(arrowstyle="fancy", color=color),
        size=size
    )

# Eggleton function ~ Roche lobe radius / orbital separation
def Egg(x):
    return (0.49*x**(2./3.))/(0.6*x**(2./3.)+np.log(1+x**(1./3.)))

f=0.8
eps_0=1.E-6

fig = plt.figure(figsize=(2*figsize,1.5*figsize))
fig1 = plt.subplot(111) # ,figsize=(1.1*figsize,figsize))

b_lab=[r'C-rich',r'O-rich']
q_lab=["",r"$q'=1$",r"$q'=10$"]
leg1=3*['']
leg1[0] = plt.Line2D((0,0),(0,1),linestyle='-',linewidth=fontsize,color=colors[0],alpha=opacities[0])
leg1[1] = plt.Line2D((0,0),(0,1),linestyle='-',linewidth=fontsize,color=colors[1],alpha=opacities[0])
# leg1[0] = plt.Line2D((0,0),(0,1),marker=mark[0],markersize=fontsize,linestyle='None',color=colors[0],markeredgewidth=1.5,markeredgecolor='k')
# leg1[1] = plt.Line2D((0,0),(0,1),marker=mark[1],markersize=fontsize,linestyle='None',color=colors[0],markeredgewidth=1.5,markeredgecolor='k')
leg2=3*['']
leg2[0] = plt.Line2D((0,0),(0,1),linestyle='-',linewidth=fontsize,color=colors[0],alpha=opacities[1])
leg2[1] = plt.Line2D((0,0),(0,1),linestyle='-',linewidth=fontsize,color=colors[1],alpha=opacities[1])
# leg2[0] = plt.Line2D((0,0),(0,1),marker=mark[0],markersize=fontsize,linestyle='None',color=colors[1],markeredgewidth=1.5,markeredgecolor='k')
# leg2[1] = plt.Line2D((0,0),(0,1),marker=mark[1],markersize=fontsize,linestyle='None',color=colors[1],markeredgewidth=1.5,markeredgecolor='k')

qq=[1.,10.]
bb=[0.1,5]
ee_min=1.4
ee_max=12.
ee=np.logspace(np.log2(ee_min),np.log2(ee_max),100,base=2.) #[0.5,0.8,1.2,2,4,8]
# ee=[4,8]

# escape speed from 1st star, in units of vorb
vesc=np.zeros(len(qq))

Rf=np.zeros((2,len(qq),len(bb),len(ee))) # 2 is for right/left
RHS=np.zeros(2)
for i in range(len(qq)):
	q=qq[i]
	vesc[i]=np.sqrt(2./((1.+1./q)*f*Egg(q)))
	print vesc[i]
	# RHS does not include centrifugal component since non-rotating star
	# in the inertial frame => initial kinetic energy in the co-rotating frame
	# cancels out w/ centrifugal component @ launch
	# Right value
	RHS[0]=(2./((1+q)*f*Egg(q)))*(q+1./(1./(f*Egg(q))-1.)) #+0.5*(1+q)*(f*Egg(q))**3.*(1./((1+q)*f*Egg(q))-1.)**2.)
	# Left value
	RHS[1]=(2./((1+q)*f*Egg(q)))*(q+1./(1./(f*Egg(q))+1.)) #+0.5*(1+q)*(f*Egg(q))**3.*(1./((1+q)*f*Egg(q))+1.)**2.)
	# print RHS[0], q, 1./(1./(f*Egg(q))-1.), 0.5*(1+q)*(f*Egg(q))**3.*(1./((1+q)*f*Egg(q))-1.)**2.
	# stop
	for j in range(len(bb)):
		b=bb[j]
		for k in range(len(ee)):
			e=ee[k]
			for p in range(len(RHS)):
				# 1st guess
				x=1.E10
				# Sanity check: 1st eps has to be > 0
				if (e**2.*(1.-1./x)**(2.*b)-RHS[p]<0.):
					# print e**2.*(1.-1./x)**(2.*b), '<', RHS[p]
					print q, b, e, '=> No free radius < 1E10'
					continue
				x1=1. # < x & (x>1)
				x2=x # > x
				while True:
	  				eps=e**2.*(1.-1./x)**(2.*b)-RHS[p]
					# print eps, x1, x, x2
					if (np.abs(eps)<eps_0):
						print q, b, e, '=>', x
						Rf[p][i][j][k]=x
						# stop
						break
					else:
						if (eps>0.):
							x2=x
							x=(x1+x)/2.
						elif (eps<0.):
							x1=x
							x=(x+x2)/2.

Rf=np.ma.masked_where(Rf<1.,Rf)
for i in range(len(qq)):
	for j in range(len(bb)):
		# print Rf[0][i][j]
		fig1.fill_between(ee/vesc[i], Rf[0][i][j], Rf[1][i][j], color=colors[j], alpha=opacities[i])

fig1.grid(which='major', linestyle='dotted', linewidth=2, alpha=1)
# fig1.grid(which='minor', linestyle='dotted', linewidth=2, alpha=1)
fig1.set_axisbelow(True)
# fig1.fill_between(np.linspace(0.,10.,100), 0, 1,color='grey')
fig1.set_xlim(0.4,8.2)
fig1.set_ylim(bottom=0.9)
# fig1.set_ylim(bottom=0.91) #,top=2.5)
fig1.set_xscale('log',basex=2)
fig1.set_yscale('log')
fig1.xaxis.set_major_formatter(FormatStrFormatter('%.1g'))
# if (ff==0.05):
# 	fig1.yaxis.set_minor_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
# elif (ff==0.8):
# fig1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: '{:g}'.format(y)))
extra = Rectangle((0, 0), 1, 1, fc="w", fill=True, edgecolor='none', linewidth=0)
legend_handle = [extra, extra, extra, extra, leg1[0], leg2[0], extra, leg1[1], leg2[1]]
leg_lab=np.concatenate([q_lab,[b_lab[0],'',''],[b_lab[1],'','']])
fig1.legend(legend_handle,leg_lab,
	loc='upper right',fontsize=fontsize,ncol=3,shadow=True,handletextpad=-1.7)
fig1.set_xlabel(r'$v_{\infty}/v_{esc}$',fontweight='bold',fontsize=fontsize)
fig1.set_ylabel(r'Free radius / Dust cond. radius', fontweight='bold', fontsize=fontsize)
fig.tight_layout()
fig.savefig(outputs+'free_radius.png',bbox_inches='tight')
# plt.show()
