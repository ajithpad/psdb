#Combined effect of BS in STN and GPe through a pcolor

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
	
#tex_lis = ['A','B','C']
tex_lis = ['D','E','F']
prefix1 = 'beta_oscil-3'
#path1 = '/data/padmanabhan/results/basal/'
path1 = '/export/data-padmanabhan/results/basal/'
prefix0 = prefix1
#log_path = '/data/padmanabhan/scripts/basal/batch/'
log_path = '/export/data-padmanabhan/scripts/basal/batch/'
log_file = open(log_path + prefix1 +'.log')
log_values = eval(log_file.read())
exc_range = log_values['exc2_ratio']
inh_range = log_values['inh2_ratio']
ext_rate = log_values['pg_rate_exc']
#add_inh_rate = log_values['extra_inh']
add_inh_rate = log_values['pg_rate_inh']
popul = 'pops_exc'
#ext_val_lis = [[4500.,2400.],[4900.,1800.],[5400.,1600.]]
ext_val_lis = [[4700.,1600.],[4900.,1800.],[5400.,2300.]]
pcol_arr = np.zeros((len(ext_val_lis),len(inh_range), len(exc_range))) 
spec_lis = []
for ext_count, ext_val in enumerate(ext_val_lis):
  for inh_count, inh_val in enumerate(inh_range):
    for exc_count, exc_val in enumerate(exc_range):
      ad1 = get_sim_3d({'pg_rate_exc':ext_val[0], 'pg_rate_inh':ext_val[1],'exc2_ratio':exc_val, 'inh2_ratio':inh_val})
      xx,xx,xx,peak_freq_gpe, peak_pow_gpe = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = popul)
      gpe_spec_val = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = popul, freq_range = [0.,1.8 * peak_freq_gpe])
      pcol_arr[ext_count, inh_count, exc_count] = gpe_spec_val

labelsize = 18. 
exc_lis = exc_range.tolist()
exc_lis.append(1.2)
exc_range = array(exc_lis)

inh_lis = inh_range.tolist()
inh_lis.append(1.2)
inh_range = array(inh_lis)
cr = arange(0.3,0.81,0.25)
fig = pl.figure(1,(12,9))
for ext_count, ext_val in enumerate(ext_val_lis):
  ax1 = fig.add_axes([0.075+(ext_count*0.28),0.1,0.225,0.35])
  pl.pcolor(exc_range, inh_range, pcol_arr[ext_count], cmap = pl.cm.coolwarm)
  pl.clim(0.3,0.8)
  
  pl.gca().set_xticks(inh_range[1::2])
  pl.gca().set_xticklabels(inh_range[1::2],size = labelsize)
  if ext_count == 1:
    pl.gca().set_xlabel('Fraction of BS neurons in STN', size = labelsize)
    
  if ext_count == 2:
    cax = fig.add_axes([0.895,0.15,0.025,0.25])
    cbar = pl.colorbar(orientation = 'vertical',cax = cax, cmap = cm.coolwarm)
    cbar.set_label('Spectral entropy of STN',position= (1.3,0.5),fontsize = labelsize, rotation = 'vertical')
    cbar.set_ticks(cr)
    pl.gca().set_yticks([])
    cticks = pl.getp(cbar.ax.axes,'yticklabels')
    pl.setp(cticks,fontsize = labelsize)
  if ext_count == 0:
    pl.gca().set_yticks([0.,0.5,1.0])
    pl.gca().set_yticklabels([0.,0.5,1.0],size = labelsize)
    pl.gca().set_ylabel('Fraction of BS neurons in GPe',size = labelsize) 
  else:
    ax1.set_yticks([])
  ax1.text(-0.1,1.2,tex_lis[ext_count], size = labelsize + 5)

#ax1 = fig.add_axes([0.7,0.1,0.225,0.35])
#for inh_count, inh_val in enumerate(inh_range):
  #spec_lis = []
  #for exc_count, exc_val in enumerate(exc_range):
      #ad1 = get_sim_3d({'pg_rate_exc':4900., 'pg_rate_inh':1800.,'exc2_ratio':exc_val, 'inh2_ratio':inh_val})
      #xx,xx,xx,peak_freq_gpe, peak_pow_gpe = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = popul)
      #gpe_spec_val = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = popul, freq_range = [0.,1.8 * peak_freq_gpe]) 
      #spec_lis.append(gpe_spec_val)
  #smooth_spec = smooth(array(spec_lis), window_len = 5)
  #ax1.plot(exc_range, smooth_spec, color = pl.cm.Greens(inh_count*25), lw = 3.)
    
  
pl.show()
  
  