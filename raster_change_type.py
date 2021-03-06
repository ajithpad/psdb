#Code to plot the rasters showing different spectral entropy

import numpy as np
import pylab as pl
import numpy as np
import itertools
import sys
sys.path.append('/users/padmanabhan/simulations/basal/scripts/common/')
import analyze_data_psd as adata
#import common.ff_analyze_data as adata
reload(adata)
import sys
from numpy import *
reload(adata)
import pdb
import pylab as pl
import matplotlib.cm as cm
import plotly.plotly as py
from plotly.graph_objs import *
import random



	
	
#prefix1 = 'beta_oscil-2'
path1 = '/media/DISK_B/results/basal/'
fil_name = 'test'

fig = pl.figure(1,(12,8))

ad1 = adata.analyze_data(fil_name,path1) 
ad1.load_spikes()
spikes = ad1.events['spikes']
labelsize = 18.

ax1 = fig.add_axes([0.075, 0.075, 0.225,0.375])
stn_spikes = spikes[spikes[:,0] < 1000.]
gpe_spikes = spikes[spikes[:,0] > 1000.]
stn_spikes_smp = random.sample(np.unique(stn_spikes[:,0]), 100)
hh = []
for ii in stn_spikes_smp:
  hh.append(stn_spikes[stn_spikes[:,0]==ii])
stn_spikes_new = np.array([item for sublist in hh for item in sublist])
gpe_spikes_smp = random.sample(np.unique(gpe_spikes[:,0]), 300)
hh = []
for ii in gpe_spikes_smp:
  hh.append(gpe_spikes[gpe_spikes[:,0]==ii])
gpe_spikes_new = np.array([item for sublist in hh for item in sublist])
#gpe_spikes_new = array(random.sample(gpe_spikes,len(gpe_spikes)/10))
ax1.plot(stn_spikes_new[:,1],stn_spikes_new[:,0], '.', color = pl.cm.Blues_r(50))
ax1.plot(gpe_spikes_new[:,1],gpe_spikes_new[:,0], '.', color = pl.cm.Greens_r(50))

psth,xx = np.histogram(stn_spikes[:,1],np.linspace(300.,900.,600))
ax2 = ax1.twinx()
ax2.plot(xx[:-1], (psth-np.mean(psth))/np.std(psth), color = pl.cm.Blues(500), lw = 3.)
ax2.plot(xx[:-1], np.zeros(len(xx[:-1])), '--', lw = 3., color = 'k', alpha = 0.7)
ax2.set_ylim(-1.25,15)
ax2.set_yticks([])

psth,xx = np.histogram(gpe_spikes[:,1],np.linspace(300.,900.,600))
ax2 = ax1.twinx()
ax2.plot(xx[:-1], (psth-np.mean(psth))/np.std(psth), color = pl.cm.Greens(500), lw = 3.)
ax2.plot(xx[:-1], np.zeros(len(xx[:-1])), '--', lw = 3., color = 'k', alpha = 0.7)
ax2.set_ylim(-7,7)
ax2.set_yticks([])

ax1.set_xlim(400,700)
ax1.set_yticks([0,1000,3500])
ax1.set_yticklabels([0,1000,3500], size = labelsize)
ax1.set_xticks(np.arange(400.,701.,300.))
ax1.set_xticklabels(np.arange(400.,701.,300.),size = labelsize)
ax1.set_ylim(0,4000)
#ax1.set_xlabel('Time(ms)', size = labelsize)
ax1.set_ylabel('Neuron id', size = labelsize, labelpad = -25)
ax1.text(590.,4000.,'D',size = labelsize + 3)
ax1.plot(ad1.pars['chg_time']*np.ones(4000), np.arange(0,4000),'--', lw = 3., color = 'y')
#ax1.text(660.,3700.,'SE = 0.3',size = labelsize, color = pl.cm.Reds_r(50))



