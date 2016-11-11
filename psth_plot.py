#Program to plot the PSTHs to show how the additional spikes change the stability of oscillations
import numpy as np
import matplotlib.pyplot as pl
import sys
sys.path.append('/users/padmanabhan/simulations/basal/scripts/common/')
import analyze_data_psd as ad
reload(ad)
results_path = '/media/DISK_B/results/basal/'
labelsize = 16
ticksize = 16
col_add = '#A6842A'
col_base =  '#246742'
col_control = '#2358B8'
fig = pl.figure(1, figsize = (8,12))
ax1 = fig.add_axes([0.11,0.35,0.4,0.2])

ad1 = ad.analyze_data('test_add_spikes',results_path)
ad1.comp_psth(form = '-', plot_fl = 1,pop_id = 'pops_inh', label = 'Add. spikes outside cycle',color = col_add)

ad1 = ad.analyze_data('test_add_spikes-control',results_path)
ad1.comp_psth(form = '-', plot_fl = 1,pop_id = 'pops_inh', label = 'Add. spikes within cycle',color = col_control)

ad1 = ad.analyze_data('test',results_path)
ad1.comp_psth(form = '-', plot_fl = 1,pop_id = 'pops_inh', label = 'No add. spikes', color = col_base)


xr = np.arange(600,761,40)
yr = np.arange(0,105,20)
#pl.legend(loc = 2, prop = {'size':14.})
pl.xlim(600,760)
pl.ylim(0,100)

#pl.gca().xaxis.set_major_locator(pl.NullLocator())
pl.gca().set_yticks(yr[1::2])
pl.gca().set_xticks(xr[::2])
pl.gca().set_yticklabels(yr[1::2], size = labelsize)
pl.gca().set_xticklabels(xr[::2], size = labelsize)
#ax1.spines['top'].set_visible(False)
#ax1.spines['bottom'].set_visible(False)
#ax1.spines['left'].set_visible(False)
#ax1.spines['right'].set_visible(False)
ax1.set_xlabel('Time (ms)',fontsize = labelsize)
ax1.set_ylabel('Firing rate (Hz)',fontsize = labelsize)
ax1.tick_params(axis = 'both', which = 'both',top = 'off', right = 'off')
#ax1.text(685,20,'D',fontsize = labelsize + 5, style = 'normal')
for tl in ax1.get_xticklabels():
    tl.set_fontsize(ticksize)
    
for tl in ax1.get_yticklabels():
    tl.set_fontsize(ticksize)



ax2 = fig.add_axes([0.57,0.35,0.4,0.2])

ad1 = ad.analyze_data('test_add_spikes-2',results_path)
ad1.comp_psth(form = '-', plot_fl = 1,pop_id = 'pops_inh', label = 'Add. spikes outside cycle',color = col_add)

ad1 = ad.analyze_data('test_add_spikes-control-2',results_path)
ad1.comp_psth(form = '-', plot_fl = 1,pop_id = 'pops_inh', label = 'Add. spikes within cycle',color = col_control)

ad1 = ad.analyze_data('test',results_path)
ad1.comp_psth(form = '-', plot_fl = 1,pop_id = 'pops_inh', label = 'No add. spikes', color = col_base)

#ad1 = ad.analyze_data('rnd_test-big_add_spikes',results_path)
#ad1.comp_psth(form = '-', plot_fl = 1, label = 'More Add. spikes outside cycle',color = np.array([0.,100.,50.])/255.)

#ad1 = ad.analyze_data('rnd_test-big_add_spikes-control',results_path)
#ad1.comp_psth(form = '-', plot_fl = 1, label = 'More Add. spikes within cycle',color = np.array([0.,50.,25.])/255.)


xr = np.arange(600,761,40)
yr = np.arange(0,105,20)
#pl.legend(loc = 2, prop = {'size':14.})
pl.xlim(600,760)
pl.ylim(0,100)
#pl.gca().xaxis.set_major_locator(pl.NullLocator())
pl.gca().set_yticks(yr[1::2])
pl.gca().set_xticks(xr[::2])
pl.gca().set_yticklabels(yr[1::2], size = labelsize)
pl.gca().set_xticklabels(xr[::2], size = labelsize)
#ax1.spines['top'].set_visible(False)
#ax1.spines['bottom'].set_visible(False)
#ax1.spines['left'].set_visible(False)
#ax1.spines['right'].set_visible(False)
ax2.set_xlabel('Time (ms)',fontsize = labelsize)
#ax1.set_ylabel('Firing rate (Hz)',fontsize = labelsize)
ax2.tick_params(axis = 'both', which = 'both',top = 'off', right = 'off')
#ax1.text(685,20,'D',fontsize = labelsize + 5, style = 'normal')
ax2.set_yticks([])
ax2.text(590,103,'D',fontsize = labelsize + 5, style = 'normal')   
ax1.text(590,103,'C',fontsize = labelsize + 5, style = 'normal')   
for tl in ax2.get_xticklabels():
    tl.set_fontsize(ticksize)

pl.legend(bbox_to_anchor = (-0.125,1.52))
fig.text(0.075, 0.975, 'A', size = labelsize + 5)
pl.show()