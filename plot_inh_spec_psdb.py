#code to plot the spectral entropy for the given values of inhibitory ratio
import numpy as np
from numpy import *
import pylab as pl
import sys
import itertools
import common.analyze_data_psd as adata
import pdb
import colorsys
from matplotlib import mpl
import matplotlib.cm as cm


prefix0 = sys.argv[1]
#prefix0 = 'mps_izh_rs_fs_ch-4'
log_path = '/media/DISK_B/scripts_backup/levcomp/batch/'
log_file = open(log_path + prefix0 +'.log')
log_values = eval(log_file.read())
exc_range = log_values['pg_rate_exc']
inh_range = log_values['extra_inh']
g_range = log_values['g']
ratio_range = log_values['inh2_ratio']
spb_range = log_values['num_spikes_burst']
#phi_range = np.array([np.round(ii,1) for ii in log_values['ff_phi']])
path1 = '/media/DISK_B/results/levcomp/'
col0 = np.array([7.,32.,32.])/255.
label_size = 13.


def change_ln(col,sc):
    ''' change lightness of color col by factor sc'''
    
    h,l,s = colorsys.rgb_to_hls(col[0],col[1],col[2])
    rgb = colorsys.hls_to_rgb(h,l*sc,s)
    #print rgb
    return array(rgb)

#col0 = change_ln(col0,1/2.)

def get_sim_3d(match_pars):
    ''' get simulation object corresponding to (g,eta,inh2_ratio)
    match_pars is a dictionary, e.g. {'g':5.6,'eta':7.6, 'inh2_ratio':0.5}
    CAUTION: check order of keys'''
    batch_pars ={'pg_rate_exc':exc_range,'num_spikes_burst':spb_range,'inh2_ratio':ratio_range,'g':g_range, 'extra_inh':inh_range}
    comb_pars = [[value for (key, value) in zip(batch_pars, values)] for values in itertools.product(*batch_pars.values())] 
    sim_id = -1
    keys_sort = batch_pars.keys() 
    match_pars_lis = [match_pars[k] for k in keys_sort]
    for ii,pars in enumerate(comb_pars):
	if match_pars_lis==np.array(pars).round(2).tolist():
	    sim_id = ii	  	    
    if sim_id>-1:
	res = adata.analyze_data(prefix0+'_%d'%sim_id,path1)
	print 'sim_id is: ',sim_id
	return res
    else:
	print 'No %s found'%match_pars.keys()
	return comb_pars,batch_pars.keys()
inh_lis = range(0,11,1)
#inh_lis = np.arange(1.,7.,1.)
col_lis = []
spb_range_new = spb_range[[0,3,5]]

def smooth(x, window_len=10, window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    import numpy as np    
    t = np.linspace(-2,2,0.1)
    x = np.sin(t)+np.random.randn(len(t))*0.1
    y = smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    
    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)
    y = np.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]

for jj,kk in enumerate(spb_range_new):
  print jj
  spec_lis = []	
  rate_lis = []
  for ii in inh_lis:
  #for ii in spb_lis:
    print ii
    ad1 = get_sim_3d({'num_spikes_burst':kk,'pg_rate_exc':352000.,'inh2_ratio':ii/10.,'g':10., 'extra_inh': 10000.})
    spec = ad1.spec_entropy(time_range = [1400.,3400.], freq_range = [0,100], pop_id = 'excitatory')
    spec_lis.append(spec)
    rate = ad1.comp_mean_rate(time_range = [1400.,3400.], pop_id = 'excitatory')
    rate_lis.append(rate)
  inh_lis_new = array(inh_lis)/10.
  col_lis.append(cm.coolwarm(jj*100))
  #Making a regression for spectral entropy
  fig1 = pl.figure(1, (8,5))
  ax1 = fig1.add_axes([0.1,0.3,0.35,0.6])
  fit1 = polyfit(inh_lis_new,spec_lis,2)
  fit1_fn = poly1d(fit1)
  #ax1.plot(inh_lis_new, fit1_fn(inh_lis_new), lw = 5., alpha = 0.8,color = change_ln(col0,jj/2.))
  smooth_spec = smooth(np.array(spec_lis), window_len = 4, window = 'hanning')
  ax1.plot(inh_lis_new, smooth_spec, '-', markersize = 3,color = cm.coolwarm(jj*100),lw = 5., alpha = 0.7, label = str(kk))
  ax1.set_xlabel('Inhibitory ratio', size = label_size)
  ax1.set_ylabel('spectral entropy', size = label_size)
  ax1.set_ylim(0.2,1.0)
  ax1.set_yticks(np.arange(0.2,1.,0.2))
  #ax2 = fig1.add_axes([0.55,0.1,0.3,0.3])
  #ax2.plot(inh_lis_new,rate_lis, lw = 5., alpha = 0.7,color = change_ln(col0,jj))
  #ax2.set_xlabel('Inhibitory ratio', size = label_size)
  #ax2.set_ylabel('Firing rate', size = label_size)
  #ax1.set_xticks(ratio_range[::2])
  #ax1.set_yticks(np.arange(0.2,1.,0.2))
  #ax2.set_xticks(ratio_range[::2])
  #ax2.set_yticks(np.arange(4.5,8.6,1.))
  
#norm = mpl.colors.Normalize(vmin = min(spb_range_new),vmax = max(spb_range_new))
#last_col = col_lis[-1]
#col_tuple = tuple(col_lis)
#cdict = {'red':[(0.0,0.0,col0[0]),
		#(1.0,last_col[0],1.0)],
	#'green':[(0.0,0.0,col0[1]),
		  #(1.0,last_col[1],1.0)],
	  #'blue': [(0.0,0.0,col0[2]),
		  #(1.0,last_col[2],1.0)]}

#my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
#mappable = pl.cm.ScalarMappable(norm = norm,cmap = my_cmap)
#mappable.set_array(spec_lis)
#cax = fig1.add_axes([0.465,0.25,0.02,0.4])
#cbar1 = fig1.colorbar(mappable,cax = cax)
#cbar1.set_ticks(spb_range_new)
#cbar1.set_label('Number of spikes per burst', size = label_size)
#for tl in ax1.get_yticklabels():
    #tl.set_fontsize(13.)
#for tl in ax1.get_xticklabels():
    #tl.set_fontsize(13.)
    
#for tl in ax2.get_yticklabels():
    #tl.set_fontsize(13.)
#for tl in ax2.get_xticklabels():
    #tl.set_fontsize(13.)
#cax = fig.add_axes([0.715,0.25,0.015,0.35])
pl.legend(loc = 'best', bbox_to_anchor =(1.2,-0.15),ncol = 3, fancybox = True)
fig1.text(0.3,0.08,'Number of spikes per burst',size = 16.)
#pl.colorbar()
pl.show()
