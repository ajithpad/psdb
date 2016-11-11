import nest
import numpy as np
import pylab as pl
import scipy.stats as st
from scipy.stats import norm
import math
import Image
import scipy

#im = Image.open('figures/sing_neu_topology.png')

#im = np.array(im).astype(np.float)/255
#height = im.size[1]
#im.resize([1,1])
tsim = 1000000.
labelsize = 15.
ticksize = 15.
#amp_arr = np.arange(64800.,65101.,100.)
#spb_arr = np.array([1.,2.,3.,4.,5.])
spb_arr = np.array([1.,3.,5.])
spk_count = []
var_count = []
burst_count = []
mean_count = []
pre_spk_count = []

def _compute_burst_index_sub(tn,limit=4,RSalpha=0.01):
    ''' compute burst_index sub procedure; called from compute_burst_index
    Computes the burst index according to the rank surprise method presented in
    Gourevitch and Eggermont 2007
    tn: spike times (ms)
    limit: smallest ISI value (ms) not to include in an identified burst
    RSalpha: significance level for discrete uniform sum distribution'''
    
    import pdb
    pdb.set_trace()
    
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
			prob = _dusd(uu,qq,N)		    
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

def _dusd(uu,qq,N):
    ''' compute discrete uniform sum distribution; required by compute_burst_index'''
    val = 0
    kmax = 1.*(uu-qq)/N
    for kk in np.arange(0,kmax+1):
	try:
	    #nom = (-1)**kk * math.factorial(uu - kk*N)
	    #denom = math.factorial(kk) * math.factorial(qq-kk) * math.factorial(uu-kk*N-qq)	    
	    #val = val + 1.*nom/denom		    
	    val = val + (-1)**kk * _cumprod2(uu-kk*N-qq+1,uu-kk*N) / (math.factorial(kk) * math.factorial(qq-kk))
	except ValueError:
	    print 'some problem'
    val = val/(N**qq)    
    return val
    
def _cumprod2(a,b):
    '''compute the cumulative product from a to b, e.g. a=2, b=5, prod2 = 2*3*4*5
    required by dusd and compute_burst_index'''	
    if a<0 or b<0:
	return 0
    else:
	return np.cumprod(np.arange(a,b+1))[-1]
	    

def compute_sing_burst_index(spikes,limit=4,RSalpha=0.01):
    '''Computes the burst index for nmax neurons and take average;
	returns burst index averaged over nmax neurons of population pop_id 
	and a detailed list of individual burst indices'''
	
    import pdb
    pdb.set_trace()
    
    burst_index_all_tot = []
    burst_index_mean_tot = []
    burst_rel_freq_tot = []
    burst_number_spikes_tot = []
    if len(spikes)==0:
	print 'Comp burst index: spike array is empty!'
	return np.nan,np.nan,np.nan,np.nan
    print 'Comp burst index'
	    

    #pop_id,spikes,neuron_nr = self.get_pop_spikes(spikes,nmax,pop_id)


    burst_index_all,burst_index_mean = _compute_burst_index_sub(spikes[:,1],limit,RSalpha)	    
    if not np.isnan(burst_index_all).any():
	burst_rel_freq = 1.*len(burst_index_all)/len(spikes)
    else:
	burst_rel_freq = np.nan
    burst_index_all_tot.append(burst_index_all)
    burst_index_mean_tot.append(burst_index_mean)
    burst_rel_freq_tot.append(burst_rel_freq)
    burst_number_spikes_tot.append((burst_rel_freq,burst_index_all,len(spikes)))
    burst_index_tot_mean = st.stats.nanmean(burst_index_mean_tot)
    burst_rel_freq_mean = st.stats.nanmean(burst_rel_freq_tot)
    tmp=[]
    for ii in burst_index_all_tot:	# flatten
	if not np.isnan(ii).any():
	    tmp.extend(ii.tolist())
    return burst_index_tot_mean*burst_rel_freq_mean,tmp,burst_rel_freq_tot,burst_number_spikes_tot
  
#fig100 = pl.figure(100)
#ax100 = fig100.add_subplot(111)
#fig101 = pl.figure(101)
#ax101 = fig101.add_subplot(111)
#col_lis = ['b','r']
#col_lis2 = ['g','k']
for count,ii in enumerate(spb_arr):
    if count > 0:
      del spikes3
    nest.ResetKernel()
    #nest.SetKernelStatus({'rng_seeds':[3]})
    volt_2 = nest.Create('voltmeter')
    dc_gen = nest.Create('poisson_generator',1,params = {'rate':70000.})
    aa = nest.Create("psdb",params = {'spb':ii})
    nest.Connect(volt_2,aa)
    nest.DivergentConnect(dc_gen,aa,weight = 1., delay = 1.)
    sd2 = nest.Create('spike_detector')
    nest.Connect(aa,sd2)
    add_poi = nest.Create('poisson_generator',params = {'rate':65000.})
    n1 = nest.Create('iaf_neuron',params = {'V_m':0.,'V_th':15000.,'E_L':0.,'V_reset':0.,'tau_m':10.})
    n2 = nest.Create('iaf_neuron',params = {'V_m':0.,'V_th':15.,'E_L':0.,'V_reset':0.,'tau_m':10.})
    vm = nest.Create('voltmeter')
    sd = nest.Create('spike_detector')
    nest.DivergentConnect(add_poi,n1+n2,weight = 1., delay = 1.)
    nest.DivergentConnect(aa,n1+n2,weight = 10., delay = 1.)
    nest.Connect(vm, n1)
    nest.Connect(n2,sd)
    nest.Simulate(tsim)
    potential = nest.GetStatus(vm, 'events')[0]
    potential2 = nest.GetStatus(volt_2, 'events')[0]
    spikes = nest.GetStatus(sd, 'events')[0]
    spikes2 = nest.GetStatus(sd2,'events')[0]
    spikes3 = np.zeros((len(spikes2['senders']),2))
    spikes3[:,1] = spikes2['times']
    spikes3[:,0] = spikes2['senders']
    #burst_ind,x,x,x = compute_sing_burst_index(spikes3)
    #burst_count.append(burst_ind)
    #ax1 = fig1.add_subplot(2,2,count+1)
    #ax1.plot(potential['times'],potential['V_m'])
    #ax2 = fig2.add_subplot(2,2,count+1)
    #ax1.plot(spikes['times'],spikes['senders'],'o')
    pre_spk_count.append((len(spikes2['senders'])/tsim)*1000)
    spk_count.append((len(spikes['senders'])/tsim)*1000)
    var_count.append(np.var(potential['V_m']))
    mean_count.append(np.mean(potential['V_m']))
    print np.mean(potential['V_m']),np.var(potential['V_m']),len(spikes['senders']), scipy.stats.shapiro(potential['V_m'])
    #ax101.plot(potential2['times'],potential2['V_m'], color = col_lis2[count])
    #ax100.hist(potential['V_m'],200, color = col_lis[count],histtype = 'step') 
    #ax101.hist(potential2['V_m'],200, color = col_lis2[count],histtype = 'step') 
fig = pl.figure(1,(8,8))
ax1 = fig.add_axes([0.09,0.075,0.8,0.375])
#ax1 = fig.add_subplot(313)
var_count_arr = np.array(var_count)
var_count_per = ((var_count_arr - var_count_arr[0])/var_count_arr[0])*100.
ax1.plot(spb_arr,var_count_per,'-o', lw = 3.,color = np.array([178.,89.,79.])/255.)
ax1.set_xlabel(' No. of spikes per burst', size = labelsize)
ax1.set_ylabel('$\%$ increase in Var. of Vm',size = labelsize,color = np.array([178.,89.,79.])/255.)
ax1.set_ylim(0.,20.)
ax1.set_yticks(np.arange(0.,21,5.))
ax1.set_xlim(1,5)
ax1.set_xticks(np.arange(1.,6.,2.))

#ax2 = fig.add_axes([0.6,0.15,0.35,0.7])
ax2 = ax1.twinx()
spk_count_arr = np.array(spk_count)
spk_count_per = ((spk_count_arr - spk_count_arr[0])/float(spk_count_arr[0]))*100.
ax2.plot(spb_arr,spk_count_per,'-o', lw = 3.,color = pl.cm.Blues(200))
ax2.set_xlabel('Number of spikes per burst',size = labelsize)
ax2.set_ylabel('$\%$ increase in firing rate',size = labelsize,color = pl.cm.Blues(200))
ax2.set_ylim(0.,50.)
ax2.set_yticks(np.arange(00.,51.,10))
ax2.set_xlim(1,5)
ax2.set_xticks(np.arange(1.,6.,2.))
ax1.text(0.7,50,'E',size = labelsize + 5)
#im = im.resize(im.size[0]/10., im.size[1]/10.)
#fig.figimage(im, 100,220 , zorder = 10)
pl.show()