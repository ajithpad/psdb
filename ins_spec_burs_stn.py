#Code to plot the GPe input against spec. entropy and burstingimport numpy as np, use the same code for STN
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
num_bars = 700
labelsize = 18.
popul = 'pops_exc'
py.sign_in('poorcricketers', 'qo8rcw0z6x')
big_spec_lis = []
fig = pl.figure(1,  (12,9))
#ax1 = fig.add_subplot(111)
#trace1 = np.zeros((len(add_inh_rate)))
data1_lis = []
data2_lis = []
#ax2 = ax1.twinx()
for ext_count, ext_val in enumerate(ext_rate[:6:2]):
  #ax1 = fig.add_axes([0.075+(ext_count*0.31),0.55,0.225,0.375])
  ax1 = fig.add_axes([0.075+(ext_count*0.31),0.55,0.225,0.35])
  ax2 = ax1.twinx()
  ax3 = ax1.twinx()
  spec_lis = []
  rate_lis = []
  gpe_rate_lis = []
  for burs_count, burs_val in enumerate(exc_range):
    ad1 = get_sim_3d({'pg_rate_exc': ext_val, 'pg_rate_inh':1700., 'exc2_ratio': burs_val, 'inh2_ratio':0.8})
    xx,xx,xx,peak_freq_gpe, peak_pow_gpe = ad1.psd(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = popul)
    gpe_freq_val = peak_freq_gpe
    rate_val = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = popul)
    rate_lis.append(rate_val)
    gpe_rate_val = ad1.comp_mean_rate(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = 'pops_inh')
    gpe_rate_lis.append(gpe_rate_val)
    gpe_spec_val = ad1.spec_entropy(time_range = [ad1.pars['T_wup'],ad1.pars['T_wup']+ ad1.pars['T_sim']],pop_id = popul, freq_range = [0.,1.8 * peak_freq_gpe])
    spec_lis.append(gpe_spec_val)
  
  smooth_gpe_rate = smooth(np.array(gpe_rate_lis), window_len = 1, window = 'hanning')
  smooth_rate = smooth(np.array(rate_lis), window_len = 1, window = 'hanning')
  smooth_spec = smooth(np.array(spec_lis), window_len = 1, window = 'hanning')
  
  #trace1  = Scatter(x = inh_range, y = smooth_spec,xaxis = 'x1', yaxis = 'y1',line = Line(color = 'rgba' + str(cm.Reds_r(ext_count * 50))))
  #data1_lis.append(trace1)
  #trace2  = Scatter(x = inh_range, y = smooth_rate,xaxis = 'x1', yaxis = 'y2',line = Line(color = 'rgba' + str(cm.Greens_r(ext_count * 50)),dash = 'dash'))
  #data1_lis.append(trace2)

  x_dum = inh_range 
  y_dum = np.arange(0.,1.1,0.1)
  z_dum = [[z]*6 for z in range(11)]
  #plt.contourf(x_dum,y_dum,z_dum, num_bars,cmap = pl.cm.Greys, alpha = 0.5)
  #bg_color = 'w'
  #plt.fill_between(x_dum,y_dum, y2 = max(y_dum), color = bg_color)
  ax1.plot(exc_range, smooth_spec, lw = 3., color = pl.cm.Reds_r(50))
  ax3.plot(exc_range, smooth_rate, lw = 3., color = pl.cm.Blues_r(50))
  ax2.plot(exc_range, smooth_gpe_rate, lw = 3., color = pl.cm.Greens_r(50))
  ax1.contourf(x_dum,y_dum,z_dum, num_bars,cmap = pl.cm.Reds, alpha = 0.5)
  ax1.set_ylim(0.,1.)
  ax2.set_ylim(26.,36.)
  ax3.set_ylim(6.,14.)
  ax1.set_xticks(np.arange(0.,1.1,0.3))
  ax1.set_xticklabels(np.arange(0.,1.1,0.3), size = labelsize)
  if ext_count == 0:
    ax1.set_yticks(np.arange(0.,1.1,0.3))
    ax1.set_yticklabels(np.arange(0.,1.1,0.3), size = labelsize, color =  pl.cm.Reds_r(50))
    ax1.set_ylabel('Spectral entropy', size = labelsize, color = pl.cm.Reds_r(50))
    
    ax2.set_yticks([])
    ax3.set_yticks([])
  if ext_count == 1:
    ax1.set_yticks([])
    ax2.set_yticks(np.arange(26.,39.,4.))
    ax2.set_yticklabels(np.arange(26.,39.,4.),size = labelsize, color = pl.cm.Greens_r(50))
    ax2.set_ylabel('Firing rate GPe(Hz)', size = labelsize, color = pl.cm.Greens_r(50))
    ax3.set_yticks([])
    ax1.set_xlabel('Frac. of BS neurons in STN', size = labelsize)
  if ext_count == 2:
    ax1.set_yticks([])
    ax3.set_yticks(np.arange(6.,15.,4.))
    ax3.set_yticklabels(np.arange(6.,15.,4.),size = labelsize, color = pl.cm.Blues_r(50))
    ax3.set_ylabel('Firing rate STN(Hz)', size = labelsize, color = pl.cm.Blues_r(50))
    ax2.set_yticks([])
  #ax1.plot(inh_range, np.ones((len(inh_range)))*0.5,'--', alpha = 0.5, color = pl.cm.Greys_r(50),lw = 3.)
  ax1.text(-0.1,1.05,tex_lis[ext_count], size = labelsize + 3)
  big_spec_lis.append(spec_lis)


#data1 = Data(data1_lis)
#layout = Layout(title = 'GPe (F = 0)',yaxis2 = YAxis(overlaying = 'y', side = 'right'))
#fig = Figure(data = data1, layout = layout)
#plot_url = py.plot(fig,filename='bursting/inp_spec_burs_stn-0')
pl.show()
#plot_url = py.plot_mpl(fig)