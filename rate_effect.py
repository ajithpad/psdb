#Code to show the rate effect on oscillations keeping the exc. and inh. ratios fixed

import numpy as np
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

textsize = 14
ticksize = 14

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

#big_gpe_ff_lis = []
#big_gpe_fr_lis = []
#big_stn_ff_lis = []
#big_stn_fr_lis = []
#big_stn_freq_lis = []
#big_stn_spec_lis = []
#big_gpe_freq_lis = []
#big_gpe_spec_lis = []

#for exc_count, exc_val in enumerate(ext_rate):
  #stn_ff_lis = []
  #stn_fr_lis = []
  #gpe_ff_lis = []
  #gpe_fr_lis = []
  #gpe_spec_lis = []
  #gpe_freq_lis = []
  #stn_spec_lis = []
  #stn_freq_lis = []
  #for inh_count, inh_val in enumerate(add_inh_rate):
    #ad1 = get_sim_3d({'pg_rate_exc':exc_val, 'pg_rate_inh': inh_val, 'exc2_ratio':0., 'inh2_ratio':0.})
    #stn_ff_val = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc', kernel_w= 2.)
    #gpe_ff_val = ad1.comp_ff(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', kernel_w= 2.)
    #stn_fr_val = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc')
    #gpe_fr_val = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
    #xx,xx,xx,peak_freq_gpe, peak_pow_gpe = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
    #gpe_freq_val = peak_freq_gpe
    #gpe_spec_val = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', freq_range = [0.,1.8 * peak_freq_gpe])
    #xx,xx,xx,peak_freq_stn, peak_pow_stn = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
    #stn_freq_val = peak_freq_stn
    #stn_spec_val = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh', freq_range = [0.,1.8 * peak_freq_stn])
	
    #stn_ff_lis.append(stn_ff_val)
    #stn_fr_lis.append(stn_fr_val)
    #gpe_ff_lis.append(gpe_ff_val)
    #gpe_fr_lis.append(gpe_fr_val)
    #gpe_spec_lis.append(gpe_spec_val)
    #gpe_freq_lis.append(gpe_freq_val)
    #stn_spec_lis.append(stn_spec_val)
    #stn_freq_lis.append(stn_freq_val)
  #big_gpe_ff_lis.append(gpe_ff_lis)
  #big_gpe_fr_lis.append(gpe_fr_lis)
  #big_stn_ff_lis.append(stn_ff_lis)
  #big_stn_fr_lis.append(stn_fr_lis)
  #big_stn_freq_lis.append(stn_freq_lis)
  #big_stn_spec_lis.append(stn_spec_lis)
  #big_gpe_freq_lis.append(gpe_freq_lis)
  #big_gpe_spec_lis.append(gpe_spec_lis)

#ext_rate.append(5600.)
#add_inh_rate.append(2600.)

#data = Data([
            #Heatmap(
                    #z = [ext_rate, add_inh_rate, big_gpe_fr_lis],colorscale = 'Earth')]) 
                    
#plot_url = py.plot(fig, filename='Earth-heatmap')

fig = pl.figure(1,(8,12))
ax1 = fig.add_axes([0.132,0.75,0.3,0.2])
pl.pcolor(add_inh_rate, ext_rate, array(big_gpe_fr_lis), cmap = pl.cm.coolwarm)
ax1.set_xlim(add_inh_rate[0], add_inh_rate[-1])
ax1.set_ylim(ext_rate[0], ext_rate[-1])
ax1.set_xticks([])
ax1.set_yticks(ext_rate[::3])
ax1.set_yticklabels(ext_rate[::3],size = textsize)
ax1.set_ylabel('Input to STN(Hz)', size = textsize)
ax1.set_title('GPe Firing rate (Hz)', size = textsize)
ax1.text(1400.,5500.,'A', size = 18)
cax = fig.add_axes([0.45,0.775,0.015,0.15])
cbar = pl.colorbar(cax = cax)
cbar.set_ticks(np.arange(25.,46.,5.))
cticks = pl.getp(cbar.ax.axes,'yticklabels')
pl.setp(cticks,fontsize = ticksize)

ax1 = fig.add_axes([0.55,0.75,0.3,0.2])
pl.pcolor(add_inh_rate, ext_rate, array(big_stn_fr_lis), cmap = pl.cm.coolwarm)
ax1.set_xlim(add_inh_rate[0], add_inh_rate[-1])
ax1.set_ylim(ext_rate[0], ext_rate[-1])
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set_title('STN Firing rate (Hz)', size = textsize)
ax1.text(1400.,5500.,'B', size = 18)
cax = fig.add_axes([0.87,0.775,0.015,0.15])
cbar = pl.colorbar(cax = cax)
cbar.set_ticks(np.arange(0.,21.,5.))
cticks = pl.getp(cbar.ax.axes,'yticklabels')
pl.setp(cticks,fontsize = ticksize)


ax1 = fig.add_axes([0.132,0.5,0.3,0.2])
pl.pcolor(add_inh_rate, ext_rate, array(big_gpe_spec_lis), cmap = pl.cm.coolwarm)
ax1.set_xlim(add_inh_rate[0], add_inh_rate[-1])
ax1.set_ylim(ext_rate[0], ext_rate[-1])
ax1.set_xticks(add_inh_rate[::4])
ax1.set_xticklabels(add_inh_rate[::4],size = textsize)
ax1.set_yticks(ext_rate[::4])
ax1.set_yticklabels(ext_rate[::4],size = textsize)
ax1.set_ylabel('Input to STN(Hz)', size = textsize)
ax1.set_xlabel('Input to GPe(Hz)', size = textsize)
ax1.set_title('GPe Spectral Entropy', size = textsize)
ax1.text(1400.,5500.,'C', size = 18)
cax = fig.add_axes([0.45,0.525,0.015,0.15])
cbar = pl.colorbar(cax = cax)
cbar.set_ticks(np.arange(0.,1.1,0.25))
cticks = pl.getp(cbar.ax.axes,'yticklabels')
pl.setp(cticks,fontsize = ticksize)


ax1 = fig.add_axes([0.55,0.5,0.3,0.2])
pl.pcolor(add_inh_rate, ext_rate, array(big_stn_spec_lis), cmap = pl.cm.coolwarm)
ax1.set_xlim(add_inh_rate[0], add_inh_rate[-1])
ax1.set_ylim(ext_rate[0], ext_rate[-1])
ax1.set_xticks(add_inh_rate[::4])
ax1.set_xticklabels(add_inh_rate[::4], size = textsize)  
ax1.set_yticks([])
ax1.set_xlabel('Input to GPe(Hz)', size = textsize)
ax1.set_title('STN Spectral Entropy', size = textsize)
ax1.text(1400.,5500.,'D', size = 18)
cax = fig.add_axes([0.87,0.525,0.015,0.15])
cbar = pl.colorbar(cax = cax)
cbar.set_ticks(np.arange(0.,1.1,0.25))
cticks = pl.getp(cbar.ax.axes,'yticklabels')
pl.setp(cticks,fontsize = ticksize)



ax1 = fig.add_axes([0.3,0.1,0.4,0.3])
pl.pcolor(array(big_gpe_fr_lis),array(big_stn_fr_lis),array(big_stn_spec_lis), cmap = cm.Reds)
ax1.set_ylabel('STN firing raate (Hz)', size = textsize)
ax1.set_xlabel('GPe firing rate (Hz)', size = textsize)
cax = fig.add_axes([0.75,0.125,0.015,0.25])
cbar = pl.colorbar(cax = cax)
cbar.set_ticks(np.arange(0.,1.1,0.25))
cticks = pl.getp(cbar.ax.axes,'yticklabels')
pl.setp(cticks,fontsize = 18)
cbar.set_label('STN spectral entropy', size = textsize)
ax1.set_xticks(np.arange(25.,46.,5.))
ax1.set_xticklabels(np.arange(25.,46.,5.),size = textsize)
ax1.set_xlim(25.,47.)
ax1.set_ylim(0.,21.)
ax1.set_yticks(np.arange(0.,21.,5.))
ax1.set_yticklabels(np.arange(0.,21.,5.),size = textsize)
ax1.text(22.,21.,'E', size = 18) 