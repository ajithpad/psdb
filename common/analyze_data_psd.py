import sys
#sys.path.append('/Users/ajith/simulations/levcomp/scripts/')
#sys.path.append('/bcfgrid/data/padmanabhan/scripts/levcomp/')
#sys.path.append('/media/DISK_B/scripts_backup/levcomp/')
sys.path.append('/Volumes/DISK_C/simulations/basal/')
import numpy as np
import pylab as pl
#import params_d as params_d
import params_d_psdb as params_d
import cPickle as cp
#import misc2
import scipy.stats as st
from scipy.stats import norm
import math
from collections import Counter
import random
import pdb


#path = '/bcfgrid/data/vlachos/results/levcomp/'

# Class with functions to analyze simulated data

class analyze_data(object):
	    
    pars = params_d.Parameters([])
    #data_path = pars.data_path
    data_path = '/Volumes/DISK_C/results/basal/' #!!!! REMEMBER to change this path if you change the results path
    #data_path = '/data/padmanabhan/results/basal/'
    def __init__(self,prefix,data_path=''):	
    	if not data_path=='':
	    self.data_path = data_path
	self.pars.prefix = prefix
	self.events = {}
	
	try:
	    self.info = self.load_info(prefix,data_path)	
	    self.pars = self.info['pars']
	    #self.info['neurons_tot'] = sum(abs(np.array(self.pars['N'])))	
	    self.pars['T_sim_range'] = [self.pars['T_wup'],self.pars['T_total']-self.pars['T_cdown']]
	except Exception,err:	    
	    raise SystemExit('\nError opening info file:\n%s !'%err)
	
    def get_pop_spikes(self,spikes,nmax,pop_id='all'):
	
	if pop_id == 'all':
	    neuron_nr = self.pars['neurons_tot']
	else:
	    #if pop_id == 'pops_exc':
	      #spiking_mode = 'spiking_mode_exc'
	    #elif pop_id == 'pops_inh':
	      #spiking_mode = 'spiking_mode_inh'
	    #pop_id = [ii for ii,string in enumerate(self.pars[spiking_mode]) if pop_id==string]
	    if pop_id==[]:
		print 'population not found'
		return [],[],[]
	    else:
		#pop_id = pop_id[0]
		if nmax==[]:
		    idx = (spikes[:,0]>=self.pars[pop_id][0]) & (spikes[:,0]<=self.pars[pop_id][-1])		    
		else: 
		    idx = (spikes[:,0]>=self.pars[pop_id][0]) & (spikes[:,0]<=self.pars[pop_id][0]+nmax)
		    neuron_nr = nmax
		neuron_nr = min(self.pars[pop_id][-1] - self.pars[pop_id][0],nmax)
		spikes = spikes[idx]	
		
	return pop_id,spikes, neuron_nr
	    
    def load_info(self,prefix,data_path=[]):
	fname = self.data_path + prefix + '.info'
	fh = open(fname,'r')
	info1 = cp.load(fh)
	return(info1)		    
	
	
    def load_spikes(self):
	if not 'spikes' in self.events.keys():
	    fname = self.data_path + self.pars['prefix'] + '_spikes.npy'
	    spikes = np.load(fname)	
	    self.events['spikes'] = spikes
	    #if len(spikes)>0:
		#print 'Imported spike data'
	    #else:
		#print 'Spike array is empty !'
	return 
		

    def load_vm(self):
	if not 'vm' in self.events.keys():
	    fname = self.data_path + self.pars['prefix'] + '_vm.npy'
	    vm = np.load(fname)	
	    self.events['vm'] = vm
	    if len(vm)>0:
		print 'Imported Vm data'
	    else:
		print 'Vm array is empty !'
	return 
	
    def plot_raster(self,nmax=100,kernel_w=1):
	''' plot raster of spikes for all neurons from 1-nr'''
	self.load_spikes()
	spikes = self.events['spikes']
	if len(spikes)==0:
	    print 'Plot raster: spike array is empty !'
	    return	
	idx = pl.find(spikes[:,0] <=nmax)
	#ax1 = pl.subplot(211)
	pl.plot(spikes[idx,1],spikes[idx,0],'.')
	#pl.xlabel('Time (ms)')
	#pl.ylabel('Neuron id')			
	#ax2 = pl.subplot(212,sharex=ax1)
	#self.comp_psth(pop_id='all',kernel_w=kernel_w)
	pl.show()    
    
    def plot_vm(self,nid=[]):
	self.load_vm()
	self.load_spikes()	
	vm = self.events['vm']
	spikes = self.events['spikes']
	if nid!=[]:
	    idx1 = vm[:,0] ==nid
	    idx2 = spikes[:,0] ==nid
	else:
	    idx1 = vm[:,0]>0
	    idx2 = idx1
	offset = np.max(vm[idx1,2])*0.8
	pl.plot(vm[idx1,1],vm[idx1,2])
	pl.plot(spikes[idx2,1],spikes[idx2,0]*0+offset,'ro')		
	pl.show()
    
    def comp_psth(self,pop_id='all',nmax=100,form = '-',kernel_w=1,res=0.1,plot_fl=1, label = 'one',kernel = 'normal', time_range = [], color=np.array([1.,1.,1.])/255.):
	''' compute psth'''
	self.load_spikes()	
	spikes = self.events['spikes']
	if len(spikes)==0:
	    print 'Comp psth: spike array is empty!'
	    return [],[]		
	
	pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	if spikes == []:
	    return [],[]	    		
	if time_range!=[]:
	    idx = (spikes[:,1]>time_range[0]) & (spikes[:,1]<=time_range[1])
	    spikes = spikes[idx]
	    sim_time = time_range[1] - time_range[0]
	
	else:
	    sim_time = self.pars['T_total']	
	
	psth,kernel,xx = misc2.kde(spikes[:,1],kernel_w,res,neuron_nr,1,kernel,sim_time)	
	
	if plot_fl:
	    pl.plot(xx,psth,form,lw = 3.,label = label, color = color)
	    pl.xlim([0,self.pars['T_sim'] + self.pars['T_wup']])
	return(psth,xx)
	
            	
    def comp_mean_rate(self,time_range=[],pop_id='all',nmax=100):
	self.load_spikes()	
	spikes = self.events['spikes']
	if len(spikes)==0:
	    print 'Comp mean rate: spike array is empty !'
	    return np.nan
	pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	if len(spikes) == 0:
	    return np.nan
	    
	if time_range!=[]:
	    idx = (spikes[:,1]>time_range[0]) & (spikes[:,1]<=time_range[1])
	    spikes = spikes[idx]
	if time_range==[]:
	  total_time = self.pars['T_total']
	else:
	  total_time = time_range[1] - time_range[0]
	mean_rate = 1.*len(spikes)/neuron_nr/total_time*1e3
	return mean_rate
	    
    def plot_rate_hist(self,bin_nr=20,pop_id = 'all',time_range = [], nmax = 100, plot_fl = 0):
	self.load_spikes()	
	spikes = self.events['spikes']
	if len(spikes)==0:
	    print 'Comp mean rate: spike array is empty !'
	    return np.nan
	pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	if len(spikes) == 0:
	    return np.nan
	    
	if time_range!=[]:
	    idx = (spikes[:,1]>time_range[0]) & (spikes[:,1]<=time_range[1])
	    spikes = spikes[idx]
	if time_range==[]:
	  total_time = self.pars['T_total']
	else:
	  total_time = time_range[1] - time_range[0]
	  
	(hh,dummy) = np.histogram(spikes[:,0],np.unique(spikes[:,0]))
	if plot_fl:
	  pl.hist(hh[:-1]/(self.pars['T_total']*1e-3),bin_nr)	# ignore last bin
	  pl.xlabel('Spike rate')
	  pl.ylabel('Counts')
	  pl.title('Spike rate histogram')
	  
	return hh,dummy
		
	
    def comp_ff(self,kernel_w=1,pop_id='all',nmax=100,time_range=[]):
	'''compute FF of population psth; use as a measure of synchrony'''	
	psth,timepoints = self.comp_psth(pop_id=pop_id,nmax=nmax,kernel_w=kernel_w,plot_fl=0)
	if len(psth)==0:
	    #print 'Comp FF: spike array is empty!'
	    return np.nan
	#print 'Comp FF'
	psth2 = psth
	if time_range:
	    idx = (timepoints>=time_range[0]) & (timepoints<=time_range[1])
	    psth2 = psth[idx]
	ff = np.var(psth2)/np.mean(psth2)
	return ff
    
    def comp_ff2(self,kernel_w=1,time_range=[]):
	'''compute FF of population psth; use as a measure of synchrony
	use spike counts per neuron, instead of population activit'''


    #def comp_osc_ind(self, pop_id='all',frange_max=1000, kernel_w=2, resolution=1e-4):
	#'''compute oscillation index'''
	
	#psth,xx1 = self.comp_psth(pop_id=pop_id,kernel_w=2,plot_fl=0)
	#if len(psth)==0:
	    #print 'Comp FF: spike array is empty!'
	    #return np.nan,np.nan
	#print 'Comp CC'
	#idx1 = (xx1 >= self.pars['T_wup']) & (xx1 <= self.pars['T_wup'] + self.pars['T_sim'])
	
	#psd,ff,vmax,fmax = misc2.psd(psth[idx1],resolution)		
	#osc_ind = vmax/sum(psd[ff<=frange_max])
    
	#return osc_ind,(psd,ff,vmax,fmax)
	
    def comp_osc_ind2(self,pop_id = 'all',nr_peaks=5,trange=5,frange_max=1000,kernel_w=2,resolution=1e-4, time_range=[]):
	'''trange (Hz): remove points +/- trange around each peak'''

	psth,xx1 = self.comp_psth(pop_id=pop_id,kernel_w=kernel_w,plot_fl=0)
	if time_range == []:
	  idx1 = (xx1 >= self.pars['T_wup']) & (xx1 <= self.pars['T_wup'] + self.pars['T_sim'])
	else:
	  idx1 = (xx1 >= time_range[0]) & (xx1 <= time_range[1])
	
	if len(psth)==0:
	    print 'Comp osc_ind: spike array is empty!'
	    return np.nan,np.nan,np.nan,np.nan

	osc_ind_lis = []
	pos_lis = []
	ff_pos1_lis = []
	ff_pos2_lis = []
	psd_lis = []
	ff_lis = []
	ff2_lis = []
	for ii in range(len(np.arange(time_range[0],time_range[1],250.))-1):
	  #import pdb
	  #pdb.set_trace()
	  vals = np.arange(time_range[0],time_range[1],250.)
	  idx1 = (xx1 >= vals[ii]) & (xx1 <= vals[ii+1])
	  if sum(idx1*1) < 10:
	    return np.nan,np.nan,np.nan,np.nan
	  else:
	    psd,ff,vmax,fmax = self.psd(time_range = time_range,pop_id = pop_id)
	    ff_step = np.diff(ff)[0]	# get frequency steps
	    trange = np.ceil(trange/ff_step)
	    psd = psd[ff<=frange_max]	# get only spectrum until frange_max
	    ff2 =ff[ff<=frange_max]
	    
	    psd=np.array(psd,dtype=float)
	    pos = np.argsort(psd)[-1::-1]   #Sorting the positions in descending order
	    psd2 = psd[pos]                 #Sorting the corresponding psds in the descending order
	    pos = np.array(pos,dtype=float)
		    
	    ii = 0
	    while ((ii-sum(np.isnan(pos[:ii]))<nr_peaks) & (ii<len(pos))):
		if not np.isnan(pos[ii]).any():	# skip nan values		
		    # find the points around peak[ii] and set to nan
		    idx = (pos>=pos[ii]-trange) & (pos<=pos[ii]+trange)
		    idx[ii] =False          # To make sure the identified peak itself is not removed
		    pos[idx] = np.nan	
		    psd2[idx] = np.nan
		ii+=1
		
	    pos = np.array(pos[np.isnan(pos).__neg__()],dtype=int) # Removing nans
	    pos = pos[:nr_peaks]	# get only nr_peaks
	    osc_ind = sum(psd[pos])/sum(psd)
	    osc_ind_lis.append(osc_ind)
	    pos_lis.append(pos)
	    if ff[pos][0] <= 200:
	      ff_pos1_lis.append(int(ff[pos][0]))
	    else:
	      ff_pos1_lis.append(0)
	    if ff[pos][1] <= 200:  
	      ff_pos2_lis.append(int(ff[pos][1]))
	    else:
	      ff_pos2_lis.append(0)
	    psd_lis.append(psd)
	    ff_lis.append(ff)
	    ff2_lis.append(ff2)
	  Counter1 = Counter(ff_pos1_lis)
	  f1,t1 = Counter1.most_common(1)[0]
	  Counter2 = Counter(ff_pos2_lis)
	  f2,t2 = Counter2.most_common(1)[0]
	  if self.pars['T_sim'] == 9000.:
	    t1 = int(t1/3)
	    t2 = int(t2/3)
	#print osc_ind,pos
	return osc_ind,psd[pos],ff[pos],psd,pos
	#return f1,f2,int(t1),int(t2),osc_ind_lis

    def comp_osc_ind_har(self,pop_id = 'all',nr_peaks=5,trange=5,frange_max=1000,kernel_w=2,resolution=1e-4, time_range=[]):
	'''trange (Hz): remove points +/- trange around each peak'''
	#To remove the harmonics

	psth,xx1 = self.comp_psth(pop_id=pop_id,kernel_w=kernel_w,plot_fl=0)
	if time_range == []:
	  idx1 = (xx1 >= self.pars['T_wup']) & (xx1 <= self.pars['T_wup'] + self.pars['T_sim'])
	else:
	  idx1 = (xx1 >= time_range[0]) & (xx1 <= time_range[1])
	
	if len(psth)==0:
	    print 'Comp osc_ind: spike array is empty!'
	    return np.nan,np.nan,np.nan,np.nan,np.nan
	
	psd,ff,vmax,fmax = self.psd(time_range = time_range,pop_id = pop_id)
	ff_step = np.diff(ff)[0]	# get frequency steps
	trange = np.ceil(trange/ff_step)
	psd = psd[ff<=frange_max]	# get only spectrum until frange_max
	ff2 =ff[ff<=frange_max]
	
	psd=np.array(psd,dtype=float)
	pos = np.argsort(psd)[-1::-1]   #Sorting the positions in descending order
	psd2 = psd[pos]                 #Sorting the corresponding psds in the descending order
	pos = np.array(pos,dtype=float)
	pos_lis_main = []	
	ii = 0
	while ((ii-sum(np.isnan(pos[:ii]))<nr_peaks) & (ii<len(pos))):
	    if not np.isnan(pos[ii]).any():	# skip nan values		
		# find the points around peak[ii] and set to nan
		idx = (pos>=pos[ii]-trange) & (pos<=pos[ii]+trange)
		idx[ii] =False          # To make sure the identified peak itself is not removed
		pos[idx] = np.nan	
		psd2[idx] = np.nan
	    ii+=1
	    
	pos = np.array(pos[np.isnan(pos).__neg__()],dtype=int) # Removing nans
	#pos = pos[:nr_peaks]	# get only nr_peaks
	osc_ind = sum(psd[pos])/sum(psd)
	pos0 = pos[0]
	pos_lis_main.append(pos0)
        for jj in pos: # Additional loop to remove harmonics and retain only the "real" peaks
	    	temp_val = float(jj)/pos0 - round(jj/pos0)
	    	if temp_val > 0 and temp_val < 0.15:
	    		pos_lis_main.append(jj)
	    	elif temp_val < 0 and abs(temp_val) > 0.85:
	    		pos_lis_main.append(jj) 
	import pdb
	pdb.set_trace()  
	print osc_ind,pos
	return osc_ind,psd[pos],ff[pos],psd[pos_lis_main],pos_lis_main
	
	
    def comp_osc_ind3(self,pop_id = 'all',nr_peaks=5,trange=5,frange_max=1000,kernel_w=2,resolution=1e-2, time_range=[], bin_w = 200.):
	'''trange (Hz): remove points +/- trange around each peak'''
	'''Identify the rate and number of spikes in each window and make a corresponding poisson process'''
	#Make the time bins
	if time_range == []:
	  time1 = 0.
	  time2 = self.pars['T_wup']+self.pars['T_sim']
	else:
	  time1 = time_range[0]
	  time2 = time_range[1]
	  
	bins = np.arange(time1,time2+(bin_w)/2,bin_w)
	osc_ind_lis = []
	pos_lis = []
	ff_pos1_lis = []
	ff_pos2_lis = []
	psd_lis = []
	ff_lis = []
	ff2_lis = []
	#Enter the bins one by one 
	for ii in (range(len(bins)-1)):
	  fir_rate = self.comp_mean_rate(pop_id = pop_id, time_range = [bins[ii],bins[ii+1]], nmax = 1000)
	  # Generate a poisson process of the same rate multiple times and make an average
	  if fir_rate > 5:    
	    Mu = fir_rate * ((bins[ii+1]-bins[ii])/1000.)
	    poiss_spikes = pl.uniform(bins[ii],bins[ii+1], int(round(Mu)))
	    poiss_spikes = poiss_spikes[poiss_spikes < bins[ii+1]]	
	    poiss_spikes.sort()
	    poiss_psth,p_kernel,xx = misc2.kde(poiss_spikes,kernel_w,resolution,1,1,'normal',bin_w) 
	    p_psd,p_ff,p_vmax,p_fmax = misc2.psd(poiss_psth, resolution)
	    p_ff_step = np.diff(p_ff)[0]
	    p_psd = p_psd[p_ff <= frange_max]
	    p_ff2 = p_ff[p_ff <= frange_max]
	
	return p_psd, p_ff2
	
	  #vals = np.arange(time1,time2,bin_w)
	  #psth,xx1 = self.comp_psth(pop_id=pop_id,kernel_w=kernel_w,plot_fl=0)
	  #idx1 = (xx1 >= vals[ii]) & (xx1 <= vals[ii+1])
	  #if sum(idx1*1) < 10:
	    #return np.nan,np.nan,np.nan,np.nan
	  #else:
	    #psd,ff,vmax,fmax = misc2.psd(psth[idx1],resolution)
	    #ff_step = np.diff(ff)[0]	# get frequency steps
	    #trange = np.ceil(trange/ff_step)
	    #psd = psd[ff<=frange_max]	# get only spectrum until frange_max
	    #ff2 =ff[ff<=frange_max]
	    
	    #psd=np.array(psd,dtype=float)
	    #pos = np.argsort(psd)[-1::-1]   #Sorting the positions in descending order
	    #psd2 = psd[pos]                 #Sorting the corresponding psds in the descending order
	    #pos = np.array(pos,dtype=float)
		    
	    #ii = 0
	    #while ((ii-sum(np.isnan(pos[:ii]))<nr_peaks) & (ii<len(pos))):
		#if not np.isnan(pos[ii]).any():	# skip nan values		
		    ## find the points around peak[ii] and set to nan
		    #idx = (pos>=pos[ii]-trange) & (pos<=pos[ii]+trange)
		    #idx[ii] =False          # To make sure the identified peak itself is not removed
		    #pos[idx] = np.nan	
		    #psd2[idx] = np.nan
		#ii+=1
		
	    #pos = np.array(pos[np.isnan(pos).__neg__()],dtype=int) # Removing nans
	    #pos = pos[:nr_peaks]	# get only nr_peaks
	    #osc_ind = sum(psd[pos])/sum(psd)
	    #osc_ind_lis.append(osc_ind)
	    #pos_lis.append(pos)
	    #if ff[pos][0] <= 200:
	      #ff_pos1_lis.append(int(ff[pos][0]))
	    #else:
	      #ff_pos1_lis.append(0)
	    #if ff[pos][1] <= 200:  
	      #ff_pos2_lis.append(int(ff[pos][1]))
	    #else:
	      #ff_pos2_lis.append(0)
	    #psd_lis.append(psd)
	    #ff_lis.append(ff)
	    #ff2_lis.append(ff2)
	  #Counter1 = Counter(ff_pos1_lis)
	  #f1,t1 = Counter1.most_common(1)[0]
	  #Counter2 = Counter(ff_pos2_lis)
	  #f2,t2 = Counter2.most_common(1)[0]
	  #if self.pars['T_sim'] == 9000.:
	    #t1 = int(t1/3)
	    #t2 = int(t2/3)
	##print osc_ind,pos
	##return osc_ind,pos,ff[pos],psd,ff,ff2
	#return f1,f2,int(t1),int(t2)
    
    def comp_cc(self,res=1,bin_size=1,time_range=[],nmax=100,pop_id = 'all', pl_flag = 0):
	''' compute pairwise correlation coefficient and the mean of all pairs'''
	self.load_spikes()
	spikes = self.events['spikes']
	
	if len(spikes)==0:
	    print 'Comp cc: spike array is empty!'
	    return np.nan,np.nan,np.nan
	print 'Comp CC'
	
	pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	
	# get only spikes within time_range
	if time_range!=[]:
	    idx = (spikes[:,1]>=time_range[0]) & (spikes[:,1]<=time_range[1])
	    spikes = spikes[idx,:]
	else:
	    time_range = [0,self.pars['T_total']]	
	    
	if len(spikes)==0:
	    print 'Comp cc: spike array is empty!'
	    return np.nan,np.nan,np.nan
	print 'Comp CC'

	ids = np.unique(spikes[:,0])	# get ids of neurons that fired    
	neuron_nr = 1
	kernel_w = 1
	tmax = self.pars['T_sim']+self.pars['T_wup']
	kernel_len = len(np.arange(-kernel_w*10,kernel_w*10+1*res,res))-1
	kde = np.zeros((nmax,tmax+kernel_len))	
	for ii,nid in enumerate(ids[:nmax]):
	    idx = spikes[:,0] == nid
	    tmp,kernel,xx = misc2.kde(spikes[idx,1],kernel_w,res,neuron_nr,1,'normal',tmax)
	    kde[ii,:] = tmp

	cc1 = np.corrcoef(kde)
	cc2 = np.triu(cc1,1).flatten()	# keep only upper triangular values	
	if pl_flag:
	  pl.hist(cc2,100,label = 'chhub')
	  pl.show()
	return cc2,st.nanmean(cc2),max(cc2)
    
	
    def comp_cv(self,pop_id='all',time_range=[],nmax=100):
	'''compute CV of ISIs for nmax neurons and take average'''
	CV_tot = []
	self.load_spikes()
	spikes = self.events['spikes']
	if len(spikes)==0:
	    print 'Comp CV: spike array is empty!'
	    return np.nan,np.nan
	print 'Comp CV'
	pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	if len(spikes) == 0:
	    print 'Comp CV: population has no spikes!'
	    return np.nan,np.nan
	if time_range!=[]:
	    idx = (spikes[:,1]>=time_range[0]) & (spikes[:,1]<=time_range[1])
	    spikes = spikes[idx,:]
	isi_tot = []
	ids = np.unique(spikes[:,0])[:nmax]
	if len(ids)==0:
	    return np.nan,np.nan
	for id_nr in ids:
	    idx = spikes[:,0] == id_nr	    
	    spikes2 = spikes[idx,1]
	    isi = np.diff(spikes2)
	    CV_tot.append(np.std(isi)/np.mean(isi))
	CV_mean = st.nanmean(CV_tot)
	return CV_mean,CV_tot
	
    def comp_cv_kl(self,pop_id='all',time_range=[],nmax=100, bin_size=0.1, plot_fl=0):
	'''compute CV of population activity (all ISIs pooled together using Kullback-Leibler divergence
	See Koyama and Shinomoto 2007'''
		
	self.load_spikes()
	spikes = self.events['spikes']
	if len(spikes)==0:
	    print 'Comp CV_KL: spike array is empty!'
	    return np.nan
	print 'Comp CV_KL'
	
	pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	if len(spikes) == 0:
	    print 'Comp CV_KL: population has no spikes!'
	    return np.nan
	if time_range!=[]:
	    idx = (spikes[:,1]>=time_range[0]) & (spikes[:,1]<=time_range[1])
	    spikes = spikes[idx,:]
	
	isi_tot = []
	ids = np.unique(spikes[:,0])[:nmax]
	if len(ids)==0:
	    return np.nan
	for id_nr in ids:
	    idx = spikes[:,0] == id_nr	    
	    spikes2 = spikes[idx,1]
	    isi = np.diff(spikes2)
	    isi_tot.extend(isi)
	if len(isi)==0:
	    return np.nan
	bins = np.arange(0.9*min(isi),1.1*max(isi),bin_size)
	hh,bins = np.histogram(isi,bins)	# normalize
	hh=1.*hh/len(isi)
	bar_w = np.mean(np.diff(bins))
		
	Hisi = -np.nansum(hh*np.log(hh))+np.log(bin_size)	#correct for discreteness
	KL = -Hisi + np.log(np.mean(isi)) + 1
	cv_kl = np.exp(-KL)
	
	if plot_fl:
	    pl.bar(bins[:-1],hh,width=bar_w)
	    
	return cv_kl
	
	
    def compute_burst_index(self,pop_id='all',nmax=100,limit=4,RSalpha=0.01): 
	'''Computes the burst index for nmax neurons and take average;
	   returns burst index averaged over nmax neurons of population pop_id 
	   and a detailed list of individual burst indices'''
	
	burst_index_all_tot = []
	burst_index_mean_tot = []
	burst_rel_freq_tot = []
	burst_number_spikes_tot = []
	self.load_spikes()
	spikes = self.events['spikes']
	if len(spikes)==0:
	    print 'Comp burst index: spike array is empty!'
	    return np.nan,np.nan,np.nan,np.nan
	print 'Comp burst index'
		

	pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	if spikes == []:
	    print 'population not found'
	    return np.nan,np.nan,np.nan,np.nan
	ids = np.unique(spikes[:,0])[:nmax]
	for id_nr in ids:
	    idx = spikes[:,0] == id_nr	    
	    spikes2 = spikes[idx,1]
	    burst_index_all,burst_index_mean = self._compute_burst_index_sub(spikes2,limit,RSalpha)	    
	    if not np.isnan(burst_index_all).any():
		burst_rel_freq = 1.*len(burst_index_all)/len(spikes2)
	    else:
		burst_rel_freq = np.nan
	    burst_index_all_tot.append(burst_index_all)
	    burst_index_mean_tot.append(burst_index_mean)
	    burst_rel_freq_tot.append(burst_rel_freq)
	    burst_number_spikes_tot.append((burst_rel_freq,burst_index_all,len(spikes2)))
	burst_index_tot_mean = st.stats.nanmean(burst_index_mean_tot)
	burst_rel_freq_mean = st.stats.nanmean(burst_rel_freq_tot)
	tmp=[]
	for ii in burst_index_all_tot:	# flatten
	    if not np.isnan(ii).any():
		tmp.extend(ii.tolist())
	return burst_index_tot_mean*burst_rel_freq_mean,tmp,burst_rel_freq_tot,burst_number_spikes_tot
	
    def _compute_burst_index_sub(self,tn,limit=4,RSalpha=0.01):
	''' compute burst_index sub procedure; called from compute_burst_index
	Computes the burst index according to the rank surprise method presented in
	Gourevitch and Eggermont 2007
	tn: spike times (ms)
	limit: smallest ISI value (ms) not to include in an identified burst
	RSalpha: significance level for discrete uniform sum distribution'''
	
	q_lim = 30	# limit for using the real distribution
	l_min = 2	# minimum length of a burst (in spikes)
	
	isi = np.diff(tn)
	N = len(isi)
	RSalpha = -np.log(RSalpha)    

	R = st.rankdata(isi)	# compute the rank of the ISIs taking into account ties

	# identify intervals of interest that may contain bursts
	isi_limit = np.diff((isi<limit)*1)
	start_int = pl.find(isi_limit==1)+1	# get index of first ISI in interval of interest
	if len(isi)==0:
	    return np.nan,np.nan
	if isi[0]<limit:
	    start_int = np.append(0,start_int)	# add first ISI if < limit
	end_int = pl.find(isi_limit==-1)
	if len(end_int)<len(start_int):
	    end_int = np.append(end_int,N)

	len_int = end_int - start_int+1		# length of intervals of interest

	#  Initializations
	archive_burst_RS = []

	# go through all intervals of interest
	for ind,nj in enumerate(start_int):
		pj = len_int[ind]
		subseq_RS = []
		# test each set of spikes
		for ii in np.arange(0,pj-(l_min-1)+1):	 
		    qq = l_min-2		# length of burst tested
		    while qq < pj-ii:
			qq = qq+1;
			uu = np.floor(sum(R[nj+ii:nj+ii+qq]))
			if qq<q_lim:	# exact discrete distribution
			    prob = self._dusd(uu,qq,N)		    
			else:
			    prob = norm.cdf((uu-qq*(N+1)/2)/np.sqrt(qq*(N**2-1)/12))
			RS = -np.log(prob)
			if RS>RSalpha:
			    subseq_RS.append((RS,ii,qq))
		if subseq_RS!=[]:
		    # sort intervals with decreasing RS value
		    subseq_RS = sorted(subseq_RS, key=lambda rsval: rsval[0], reverse=True)
		    tmp = np.array(subseq_RS)
		    while len(tmp)>0:
			current_burst = tmp[0]	# extract most surprising burst
			archive_burst_RS.append((current_burst[0],current_burst[2]+1,nj+current_burst[1]))		
			idx = np.array(1*(tmp[:,1]+tmp[:,2]-1<current_burst[1]) | 1*(tmp[:,1]>current_burst[1]+current_burst[2]-1),dtype=bool)		
			tmp = tmp[idx,:]
	    
	tmp = np.array(archive_burst_RS)  
	if len(tmp):
	    idx = tmp[:,0]!=np.inf	# remove inf values
	    tmp = tmp[idx]	
	    return tmp[:,0],np.mean(tmp[:,0])
	else:
	    return np.array(np.nan),np.array(np.nan)
    
    def _dusd(self,uu,qq,N):
	''' compute discrete uniform sum distribution; required by compute_burst_index'''
	val = 0
	kmax = 1.*(uu-qq)/N
	for kk in np.arange(0,kmax+1):
	    try:
		#nom = (-1)**kk * math.factorial(uu - kk*N)
		#denom = math.factorial(kk) * math.factorial(qq-kk) * math.factorial(uu-kk*N-qq)	    
		#val = val + 1.*nom/denom		    
		val = val + (-1)**kk * self._cumprod2(uu-kk*N-qq+1,uu-kk*N) / (math.factorial(kk) * math.factorial(qq-kk))
	    except ValueError:
		print 'some problem'
	val = val/(N**qq)    
	return val
	
    def _cumprod2(self,a,b):
	'''compute the cumulative product from a to b, e.g. a=2, b=5, prod2 = 2*3*4*5
	required by dusd and compute_burst_index'''	
	if a<0 or b<0:
	    return 0
	else:
	    return np.cumprod(np.arange(a,b+1))[-1]
	    
    
    def comp_sta(self, pop_type = 'fast_spiking',resolution = 0.1,time_min = 200, time_max = 200, neuron_nr = 4210.):
	
	self.load_spikes()
	spikes = self.events['spikes']
	psth,xx = self.comp_psth(pop_id = 'regular_spiking',res = resolution,plot_fl=0)
	
	dummy,trigger_spikes,dummy = self.get_pop_spikes(spikes,nmax = [],pop_id = pop_type)
	if (trigger_spikes == []):
	  return np.nan
	else:
	  z = trigger_spikes[:,1][(trigger_spikes[:,0]== neuron_nr)]
	  if len(z) == 0:
	    return np.nan
	  else:
	    sta_list = []
	    for i in range(0,len(z)):
	      x= z[i]+time_max
	      y= z[i]-time_min
	      
	      d = psth[(xx >= y) & (xx <= x)]
	      print len(d)
	      if len(d)==4000:
		sta_list.append(d)
	      
	    fin_st_sum = sum(sta_list)
	    fin_sta = fin_st_sum/len(z)
	    tt = np.arange(-time_min,time_max,resolution)	    
	return fin_sta,tt
	
	
	
    def compute_sing_burst_index(self,nmax=100,limit=4,RSalpha=0.01,nid = 4021):
	'''Computes the burst index for nmax neurons and take average;
	   returns burst index averaged over nmax neurons of population pop_id 
	   and a detailed list of individual burst indices'''
	
	burst_index_all_tot = []
	burst_index_mean_tot = []
	burst_rel_freq_tot = []
	burst_number_spikes_tot = []
	self.load_spikes()
	spikes = self.events['spikes']
	spikes = spikes[spikes[:,0]==nid]
	if len(spikes)==0:
	    print 'Comp burst index: spike array is empty!'
	    return np.nan,np.nan,np.nan,np.nan
	print 'Comp burst index'
		

	#pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	if spikes == []:
	    print 'population not found'
	    return np.nan,np.nan,np.nan,np.nan
	ids = np.unique(spikes[:,0])
	for id_nr in ids:
	    idx = spikes[:,0] == id_nr	    
	    spikes2 = spikes[idx,1]
	    burst_index_all,burst_index_mean = self._compute_burst_index_sub(spikes2,limit,RSalpha)	    
	    if not np.isnan(burst_index_all).any():
		burst_rel_freq = 1.*len(burst_index_all)/len(spikes2)
	    else:
		burst_rel_freq = np.nan
	    burst_index_all_tot.append(burst_index_all)
	    burst_index_mean_tot.append(burst_index_mean)
	    burst_rel_freq_tot.append(burst_rel_freq)
	    burst_number_spikes_tot.append((burst_rel_freq,burst_index_all,len(spikes2)))
	burst_index_tot_mean = st.stats.nanmean(burst_index_mean_tot)
	burst_rel_freq_mean = st.stats.nanmean(burst_rel_freq_tot)
	tmp=[]
	for ii in burst_index_all_tot:	# flatten
	    if not np.isnan(ii).any():
		tmp.extend(ii.tolist())
	return burst_index_tot_mean*burst_rel_freq_mean,tmp,burst_rel_freq_tot,burst_number_spikes_tot
	
	
    	
    def compFanoFactor(self,time_range=[],pop_id= 'all',nmax=100):
	''' Compute the fano factor for the spike trains given in sp'''
	self.load_spikes()	
	spikes = self.events['spikes']
	if len(spikes)==0:
	    print 'Comp mean rate: spike array is empty !'
	    return np.nan
	pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	if len(spikes) == 0:
	    return np.nan
	    
	if time_range!=[]:
	    idx = (spikes[:,1]>time_range[0]) & (spikes[:,1]<=time_range[1])
	    spikes = spikes[idx]
	if time_range==[]:
	  total_time = self.pars['T_total']
	else:
	  total_time = time_range[1] - time_range[0]
	ids = np.unique(spikes[:,0])[:nmax]

	counts = np.zeros((len(ids),))
	for i in np.arange(len(ids)):
		counts[i] = len(spikes[spikes[:,0]==i,:])

	FF = (st.nanstd(counts))**2/st.nanmean(counts)
	return FF
	
    def coll_pg_rate(self):
	pg = self.info['pars']['pg_rate']
	return pg
	

	
    def comp_ind_rate(self,pop_id = 'all',time_range=[],nmax=100):
      '''compute rate of individual neurons for nmax neurons and arrange them in order and take average'''
      Rate_tot = []
      self.load_spikes()
      spikes = self.events['spikes']
      if len(spikes)==0:
	  print 'Comp Rate: spike array is empty!'
	  return np.nan,np.nan
      print 'Comp Rate'
      pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
      if len(spikes) == 0:
	  print 'Comp Rate: population has no spikes!'
	  return np.nan,np.nan
      if time_range!=[]:
	  idx = (spikes[:,1]>=time_range[0]) & (spikes[:,1]<=time_range[1])
	  spikes = spikes[idx,:]
      
      if time_range==[]:
	  total_time = self.pars['T_total']
      else:
	total_time = time_range[1] - time_range[0]
 
      ids = np.unique(spikes[:,0])[:nmax]
      if len(ids)==0:
	  return np.nan,np.nan
      dtype = [('id',float),('rate',float)]
      for id_nr in ids:
	  idx = spikes[:,0] == id_nr	    
	  spikes2 = spikes[idx,1]
	  rate_ind = (len(spikes2)/total_time)*1000
	  rate_ind = (id_nr,rate_ind)
	  Rate_tot.append(rate_ind)
      Rate_tot = np.array(Rate_tot, dtype = dtype)
      sort_rate = np.sort(Rate_tot, order = 'rate')[::-1]
      sort_rate = sort_rate[:int(neuron_nr/10)]
      Rate_mean = st.nanmean(sort_rate['rate'])
      return Rate_mean,sort_rate
      
    def comp_fil_ff(self,pop_id = 'all',time_range=[], nmax = 500, kernel_w = 2., res = 0.1):
      
      Rate_mean, sort_rate = self.comp_ind_rate(pop_id = pop_id, time_range = time_range, nmax = nmax)
      neuron_id = sort_rate['id']
      self.load_spikes()
      spikes = self.events['spikes']
      spike_lis = []
      dtype = [('id',float),('time',float)]
      for ii in neuron_id:
	idx = (spikes[:,0]==ii)
	spike_neu = spikes[idx]
	spike_lis.append(np.array(spike_neu))
      sort_spike = np.concatenate(spike_lis)
      neuron_nr = len(np.unique(sort_spike[:,0]))
      spikes = sort_spike
      if spikes == []:
	  return [],[]	
      if time_range!= []:
	  idx = (spikes[:,1]>time_range[0]) & (spikes[:,1]<time_range[1])
	  spikes = spikes[idx,:]
	  sim_time = time_range[1] - time_range[0]	
      elif time_range == []:
	  sim_time = self.pars['T_total']
      psth,kernel,xx = misc2.kde(spikes[:,1],kernel_w,res,neuron_nr,1,'normal',sim_time)
      if len(psth)==0:
	  print 'Comp FF: spike array is empty!'
	  return np.nan
      print 'Comp FF'
      ff = np.var(psth)/np.mean(psth)
      return sort_spike,spike_lis,ff
      
    # To compute the power spectral density and return the two biggest values 
    def psd(self, bin_w = 5., nmax = 4000, time_range = [], pop_id = 'all'):
      self.load_spikes()
      spikes = self.events['spikes']
      pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
      if len(spikes) == 0:
	print 'psd: spike array is empty'
	return np.nan, np.nan, np.nan, np.nan
      if time_range!=[]:
	    idx = (spikes[:,1]>time_range[0]) & (spikes[:,1]<=time_range[1])
	    spikes = spikes[idx]
      if time_range==[]:
	total_time = self.pars['T_total']
      else:
	total_time = time_range[1] - time_range[0]
      
      if len(spikes) == 0:
	print 'psd: spike array is empty'
	return np.nan, np.nan, np.nan, np.nan
	
      ids = np.unique(spikes[:,0])[:nmax]
      nr_neurons = len(ids)
      #psd, max_value, freq,h = misc2.psd_sp(spikes[:,1],nr_bins,nr_neurons)
      bins = np.arange(time_range[0],time_range[1],bin_w)
      a,b = np.histogram(spikes[:,1], bins)
      ff = abs(np.fft.fft(a- np.mean(a)))**2
      Fs = 1./(bin_w*0.001)
      freq2 = np.fft.fftfreq(len(bins))[0:len(bins/2)+1]
      freq = np.linspace(0,Fs/2,len(ff)/2+1)
      px = ff[0:len(ff)/2+1]
      max_px = np.max(px[1:])
      idx = px == max_px
      corr_freq = freq[pl.find(idx)]
      new_px = px
      max_pow = new_px[pl.find(idx)]
      return new_px,freq, freq2, corr_freq[0], max_pow
      
    def inh_mean_rate(self,time_range = [], nmax = 100):
	fs_rate = self.comp_mean_rate(time_range = time_range, pop_id = 'fast_spiking', nmax = nmax)
	ch_rate = self.comp_mean_rate(time_range = time_range, pop_id = 'chattering', nmax = nmax)
	if np.isnan(fs_rate):
	  fs_rate = 0.
	if np.isnan(ch_rate):
	  ch_rate = 0.
	inh_ratio = self.pars['inh2_ratio']
	mean_rate = (fs_rate * (1-inh_ratio)) + (ch_rate * inh_ratio)
	return mean_rate 
	
	
    def compute_sing_burst_index_one_neuron(self,limit=4,RSalpha=0.01,nid = 1,spikes = []):
	'''Computes the burst index for nmax neurons and take average;
	   returns burst index averaged over nmax neurons of population pop_id 
	   and a detailed list of individual burst indices'''
	
	burst_index_all_tot = []
	burst_index_mean_tot = []
	burst_rel_freq_tot = []
	burst_number_spikes_tot = []
	if len(spikes)==0:
	    print 'Comp burst index: spike array is empty!'
	    return np.nan,np.nan,np.nan,np.nan
	print 'Comp burst index'
		

	#pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)
	if spikes == []:
	    print 'population not found'
	    return np.nan,np.nan,np.nan,np.nan	    
	spikes2 = spikes[:,1]
	burst_index_all,burst_index_mean = self._compute_burst_index_sub(spikes2,limit,RSalpha)	    
	if not np.isnan(burst_index_all).any():
	    burst_rel_freq = 1.*len(burst_index_all)/len(spikes2)
	else:
	    burst_rel_freq = np.nan
	burst_index_all_tot.append(burst_index_all)
	burst_index_mean_tot.append(burst_index_mean)
	burst_rel_freq_tot.append(burst_rel_freq)
	burst_number_spikes_tot.append((burst_rel_freq,burst_index_all,len(spikes2)))
	burst_index_tot_mean = st.stats.nanmean(burst_index_mean_tot)
	burst_rel_freq_mean = st.stats.nanmean(burst_rel_freq_tot)
	tmp=[]
	for ii in burst_index_all_tot:	# flatten
	    if not np.isnan(ii).any():
		tmp.extend(ii.tolist())
	return burst_index_tot_mean*burst_rel_freq_mean,tmp,burst_rel_freq_tot,burst_number_spikes_tot

	
    def spec_entropy(self,bin_w = 5.,pop_id = 'pops_exc',time_range=[],freq_range = []):
      '''Function to calculate the spectral entropy'''
      power,freq,dummy,dummy,dummy = self.psd(pop_id = pop_id,bin_w = bin_w,time_range = time_range)
      if freq_range != []:
	power = power[(freq>freq_range[0]) & (freq < freq_range[1])]
	freq = freq[(freq>freq_range[0]) & (freq < freq_range[1])]	
      k = len(freq)
      power = power/sum(power)
      sum_power = 0
      for ii in range(k):
	sum_power += (power[ii]*np.log(power[ii]))
      spec_ent = -(sum_power/np.log(k))  
      return spec_ent
      
    def calc_max_psth(self,pop_id = 'regular_spiking', time_range = [1020.,1060.]):
      '''Function to calculate the height of the peak in a given PSTH time range'''
      
      psth = self.comp_psth(pop_id = pop_id, plot_fl = 0)
      low_range = time_range[0]
      high_range = time_range[1]
      low_psth1 = psth[1][psth[1] > low_range]
      low_psth0 = psth[0][psth[1] > low_range]
      fin_psth1 = low_psth1[low_psth1 < high_range]
      fin_psth0 = low_psth0[low_psth1 < high_range]
      
      norm_max = max(fin_psth0)/sum(fin_psth0)
      return fin_psth1, fin_psth0, max(fin_psth0),norm_max 
      
    
    def find_peaks(self,pop_id = 'regular_spiking',thr=15.,time_range = []):
      #To find the peaks in a given time range from a PSTH
      
      trace,timepoints = self.comp_psth(pop_id = 'regular_spiking',plot_fl = 0)
      peak_time, peak_value = misc2.find_peaks(timepoints,trace,thr,time_range)
      
      return peak_time[0],peak_value[0]      
      
      
	  
	



	
            

	
	    
