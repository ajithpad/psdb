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


p1 = []
lab_names = []
legend_dict = dict()
def get_sim_3d(match_pars):
    ''' get simulation object corresponding to (g,eta,inh2_ratio)
    match_pars is a dictionary, e.g. {'g':5.6,'eta':7.6, 'inh2_ratio':0.5}
    CAUTION: check order of keys'''  
    #batch_pars ={'pg_rate_exc':ext_rate,'extra_inh':add_inh_rate,'exc2_ratio':exc_range,'inh2_ratio':inh_range }
    batch_pars ={'pg_rate_exc':ext_rate,'pg_rate_inh':add_inh_rate,'exc2_ratio':exc_range,'inh2_ratio':inh_range }
    comb_pars = [[value for (key, value) in zip(batch_pars, values)] for values in itertools.product(*batch_pars.values())] 
    sim_id = -1
    keys_sort = batch_pars.keys() 
    match_pars_lis = [match_pars[k] for k in keys_sort]
    for ii,pars in enumerate(comb_pars):
	if match_pars_lis==np.array(pars).round(2).tolist():
	    sim_id = ii	  	
	    #print pars
    if sim_id>-1:
	res = adata.analyze_data(prefix0+'_%d'%sim_id,path1)
	#print 'sim_id is: ',sim_id
	return res
    else:
	print 'No %s found'%match_pars.keys()
	return comb_pars,batch_pars.keys()
	
	
prefix1 = 'beta_oscil-2'
path1 = '/data/padmanabhan/results/basal/'
prefix0 = prefix1
log_path = '/data/padmanabhan/scripts/basal/batch/'
log_file = open(log_path + prefix1 +'.log')
log_values = eval(log_file.read())
exc_range = log_values['exc2_ratio']
inh_range = log_values['inh2_ratio']
ext_rate = log_values['pg_rate_exc']
#add_inh_rate = log_values['extra_inh']
add_inh_rate = log_values['pg_rate_inh']

fig = pl.figure(1,(12,8))
#For SE < 0.3
#ad1 = get_sim_3d({'pg_rate_exc':5200.,'pg_rate_inh':1700.,'inh2_ratio':0.,'exc2_ratio':0.})
#fil_name = 'se_3_spikes.npy'
#results_path = '/media/DISK_B/results/basal/'
ad1 = get_sim_3d({'pg_rate_exc':5200.,'pg_rate_inh':1700.,'inh2_ratio':0.,'exc2_ratio':0.})
ad1.load_spikes()
spikes = ad1.events['spikes']
labelsize = 18.
#spikes = np.load(results_path + fil_name)
#ad1.load_spikes()
#spikes = ad1.events['spikes']
ax1 = fig.add_axes([0.075, 0.075, 0.225,0.375])
stn_spikes = spikes[spikes[:,0] < 1000.]
gpe_spikes = spikes[spikes[:,0] > 1000.]
ax1.plot(stn_spikes[:,1],stn_spikes[:,0], '.', color = pl.cm.Blues_r(50))
ax1.plot(gpe_spikes[:,1],gpe_spikes[:,0], '.', color = pl.cm.Greens_r(50))

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

ax1.set_xlim(600,700)
ax1.set_yticks([0,1000,3500])
ax1.set_yticklabels([0,1000,3500], size = labelsize)
ax1.set_xticks(np.arange(600.,701.,100.))
ax1.set_xticklabels(np.arange(600.,701.,100.),size = labelsize)
ax1.set_ylim(0,4000)
#ax1.set_xlabel('Time(ms)', size = labelsize)
ax1.set_ylabel('Neuron id', size = labelsize, labelpad = -25)
ax1.text(590.,4000.,'D',size = labelsize + 3)
ax1.text(660.,3700.,'SE = 0.3',size = labelsize, color = pl.cm.Reds_r(50))



ad1 = get_sim_3d({'pg_rate_exc':5200.,'pg_rate_inh':1900.,'inh2_ratio':0.,'exc2_ratio':0.})
ad1.load_spikes()
spikes = ad1.events['spikes']
#labelsize = 18.
#spikes = np.load(results_path + fil_name)
#ad1.load_spikes()
#spikes = ad1.events['spikes']
ax1 = fig.add_axes([0.385, 0.075, 0.225,0.375])
stn_spikes = spikes[spikes[:,0] < 1000.]
gpe_spikes = spikes[spikes[:,0] > 1000.]
ax1.plot(stn_spikes[:,1],stn_spikes[:,0], '.', color = pl.cm.Blues_r(50))
ax1.plot(gpe_spikes[:,1],gpe_spikes[:,0], '.', color = pl.cm.Greens_r(50))

psth,xx = np.histogram(stn_spikes[:,1],np.linspace(600.,1200.,600))
ax2 = ax1.twinx()
ax2.plot(xx[:-1], (psth-np.mean(psth))/np.std(psth), color = pl.cm.Blues(500), lw = 3.)
ax2.plot(xx[:-1], np.zeros(len(xx[:-1])), '--', lw = 3., color = 'k', alpha = 0.7)
ax2.set_ylim(-1.25,15)
ax2.set_yticks([])

psth,xx = np.histogram(gpe_spikes[:,1],np.linspace(600.,1200.,600))
ax2 = ax1.twinx()
ax2.plot(xx[:-1], (psth-np.mean(psth))/np.std(psth), color = pl.cm.Greens(500), lw = 3.)
ax2.plot(xx[:-1], np.zeros(len(xx[:-1])), '--', lw = 3., color = 'k', alpha = 0.7)
ax2.set_ylim(-7,7)
ax2.set_yticks([])
ax1.set_xlim(1000,1100)
ax1.set_yticks([])
#ax1.set_yticklabels([], size = labelsize)
ax1.set_xticks(np.arange(1000.,1101.,100.))
ax1.set_xticklabels(np.arange(600.,701.,100.),size = labelsize)
ax1.set_ylim(0,4000)
ax1.set_xlabel('Time(ms)', size = labelsize)
ax1.text(990.,4000.,'E',size = labelsize + 3)
ax1.text(1060.,3700.,'SE = 0.5',size = labelsize, color = pl.cm.Reds_r(50))

ad1 = get_sim_3d({'pg_rate_exc':5200.,'pg_rate_inh':2200.,'inh2_ratio':0.,'exc2_ratio':0.6})
ad1.load_spikes()
spikes = ad1.events['spikes']
#labelsize = 18.
#spikes = np.load(results_path + fil_name)
#ad1.load_spikes()
#spikes = ad1.events['spikes']
ax1 = fig.add_axes([0.695, 0.075, 0.225,0.375])
stn_spikes = spikes[spikes[:,0] < 1000.]
gpe_spikes = spikes[spikes[:,0] > 1000.]
ax1.plot(stn_spikes[:,1],stn_spikes[:,0], '.', color = pl.cm.Blues_r(50))
ax1.plot(gpe_spikes[:,1],gpe_spikes[:,0], '.', color = pl.cm.Greens_r(50))
ax1.set_xlim(1000,1100)
ax1.set_yticks([])
psth,xx = np.histogram(stn_spikes[:,1],np.linspace(600.,1200.,600))
ax2 = ax1.twinx()
ax2.plot(xx[:-1], (psth-np.mean(psth))/np.std(psth), color = pl.cm.Blues(500), lw = 3.)
ax2.plot(xx[:-1], np.zeros(len(xx[:-1])), '--', lw = 3., color = 'k', alpha = 0.7)
ax2.set_ylim(-1.25,15)
ax2.set_yticks([])

psth,xx = np.histogram(gpe_spikes[:,1],np.linspace(600.,1200.,600))
ax2 = ax1.twinx()
ax2.plot(xx[:-1], (psth-np.mean(psth))/np.std(psth), color = pl.cm.Greens(500), lw = 3.)
ax2.plot(xx[:-1], np.zeros(len(xx[:-1])), '--', lw = 3., color = 'k', alpha = 0.7)
ax2.set_ylim(-7,7)
ax2.set_yticks([])
ax1.set_xlim(1000,1100)
ax1.set_yticks([])

#ax1.set_yticklabels([], size = labelsize)
ax1.set_xticks(np.arange(1000.,1101.,100.))
ax1.set_xticklabels(np.arange(600.,701.,100.),size = labelsize)
ax1.set_ylim(0,4000)
#ax1.set_xlabel('Time(ms)', size = labelsize)
ax1.text(990.,4000.,'F',size = labelsize + 3)
ax1.text(1060.,3700.,'SE = 0.8',size = labelsize, color = pl.cm.Reds_r(50))