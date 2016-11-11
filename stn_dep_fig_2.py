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
	
tex_lis = ['A','B','C']
#tex_lis = ['D','E','F']
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
popul = 'pops_exc'
#ext_val_lis = [[4500.,2400.],[5100.,1900.],[5400.,1600.]]
#pcol_arr = np.zeros((len(ext_val_lis),len(inh_range), len(exc_range))) 
spec_lis = []
stn_fr_lis = []
gpe_fr_lis = []

for ii in range(len(ext_rate)):
  print ext_rate[len(ext_rate)-2-ii],add_inh_rate[ii] 
  ad1 = get_sim_3d({'pg_rate_exc':ext_rate[len(ext_rate)-2-ii],'pg_rate_inh':add_inh_rate[ii],'exc2_ratio':0., 'inh2_ratio':0.})
  gpe_fr = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
  stn_fr = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc')
  xx,xx,xx,peak_freq_stn, peak_pow_stn = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_exc')
  stn_freq_val = peak_freq_stn
  stn_spec_val = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = popul, freq_range = [0.,1.8 * peak_freq_stn])
  spec_lis.append(stn_spec_val)
  gpe_fr_lis.append(gpe_fr)
  stn_fr_lis.append(stn_fr)



