#! /usr/bin/env python

from pylab import *
import shutil
import shelve
#import tables
import nest
import numpy
import numpy.random as nprandom
import pdb
#import NeuroTools	 
#from NeuroTools.parameters import ParameterSet
#from NeuroTools.parameters import ParameterRange
import time
import os
#import sqlite3
#import mako
#from mako.template import Template
#import webbrowser

def pb_stdp (p, max_time, sim_res):

	tau = p.tau_plus			#assuming that tau_plus = tau_minus

	pre = load('spikes1.dat')
	post = load('spikes2.dat')

	#ts = arange(0,20,3)
	f = lambda    t :  exp(-t/tau)							#  define exponential kernel
	g = lambda t,ti :  f(t-ti) * ones(size(t))*(t>=ti)		# place exponential kernel at the occurence of spikes

	xx = arange(0,max_time,0.1)
	val1 = empty([size(pre),size(xx)])
	val2 = empty([size(post),size(xx)])

	figure()


	for i in arange(size(pre)):			# apply an exponential kernel at the position of each presynaptic spike

		val1[i,:] = g(xx,pre[i])

	for i in arange(size(post)):			# apply an exponential kernel at the position of each postsynaptic spike

		val2[i,:] = g(xx,post[i])	

	pre_trace = cumsum(val1,0)[-1]			# add all kernels to get presynaptic trace (see Morrison et al. 2008, fig.3) ;  get last row 
	post_trace = cumsum(val2,0)[-1]			# add all kernels to get postsynaptic trace (see Morrison et al. 2008, fig.3) ;  get last row 


	weights1 = zeros(max_time/sim_res)
	weights2 = zeros(max_time/sim_res)

	for i in pre:						#  negative changes at presynaptic spike times proportional to post_trace and weight function F-
		
		idx = i/sim_res
		weights1[idx : ] = weights1[idx : ] + p.alpha * p.A_plus * post_trace[idx] * p.w_max

			
	for i in post:						#  positive changes at postsynaptic spike times proportional to pre_trace and weight function F+
		
		idx = i/sim_res
		weights2[idx : ] = weights2[idx : ] + p.A_plus * pre_trace[idx] * p.w_max
		
	weights = weights2 -  weights1
	weights = weights + p.pop1_pop2_w	# add initial weight to see fluctuations around it

		
	subplot(3,1,1)	
	stem(pre,ones(size(pre)),'r',markerfmt='w.')
	plot(xx,pre_trace)

	subplot(3,1,2)	
	stem(post,ones(size(post)),'b',markerfmt='w.')
	plot(xx,post_trace,'r')


	subplot(3,1,3)
	plot(xx,weights,'k')
	# plot(xx,weights1,'r')
	#plot(xx,weights2,'b')
	x1,x2 = xlim()
	hlines(p.pop1_pop2_w,x1,x2)
	ylim(ylim()[0] * 0.9, ylim()[1]*1.1)
	show()
	

def savefigure():
	nfig = load('/home/vlachos/model/figures/log')
	nfig += 1
	save('/home/vlachos/model/figures/figurecount.txt',nfig)
	fname = '/home/vlachos/model/figures/stdp2-0'+str(int(nfig))+'.png'
	savefig(fname)
	source = '/home/vlachos/model/param02.py'
	target = '/home/vlachos/model/figures/stdp2-0'+str(int(nfig))+'.par'
	shutil.copy(source,target)

def common_elements(v1,v2):
	"""Returns a vector of size v1 that is true for any value in v1 that matches any value in v2."""
	# v1 and v2 have to be arrays
	
	common = zeros(len(v1),'bool')
	for i in v2:
		common[v1==i]=True
	
	return common
	
def activity_matrix(v1,v2,dims):
	""" Returns a binary valued matrix with dimensions dims. The elements of the matrix have a value of 
    one if the corresponding elements in v1  appear also in v2. The rest is zero."""
	
	offset = v1[0]
	v1 = v1 - offset
	v2 = v2 - offset
	ind = common_elements(v1,v2)
	v1 = v1*0
	v1[ind] = 1
	v1 = v1.reshape(dims)
	return v1

def slw(vlen, window, overlap):
	
	"""Returns the edges for a sliding window to be applied in a vector of length vlen. 
	window is the window size and overlap the percentage of overlaping windows.
	If a window exceeds the range then it's truncated"""
	
	wstep = round((1-overlap) * window)     # the step size for each window
	idx = arange(0,vlen-window/2,wstep)     # the points around which each window is centered
	edges = empty([len(idx),2])       
	for i in arange(len(idx)):
	    edges[i,:] = [idx[i],min(vlen,idx[i]+window-1)] # truncate window if it exceeds the range
	
	return edges
	
def get_idx(v1,edges):
	"""Returns a list that has len(edges) arrays as elements. The values of each array are true for every
	v1 than falls within edges[i], otherwise it is false."""
	
	idx = []
	for i in arange(len(edges)):
		idx = idx + [(v1 >= edges[i,0]) & (v1 < edges[i,1])]
		
	return idx

#def mean_v(mat):					
#	'''Computes voltage per unique time point averaged over all neurons
#	mat has 3 columns, 0: id, 1: time-point, 2:voltage-value 
#	'''
#	tp = unique(mat[:,1])			# get unique time points
#	val = empty((len(tp),2))		# allocate memory for each time point
#	rows = 0
#	for i in tp:
#		idx = mat[:,1] == i
#		val[rows,0] = i
#		val[rows,1] = mean(mat[idx,2])	# calc average over all neurons for each unique time point	
#		rows +=1
#	return val
	
def mean_v(mat):					
	'''Computes voltage per unique time point averaged over all neurons.
	mat has 3 columns, 0: id, 1: time-point, 2:voltage-value 
	'''
	tp = unique(mat[:,1])			# get unique time points
	ids = unique(mat[:,0])			# get unique ids
	l = size(ids)	
	val = empty((len(tp),2))		# allocate memory for each time point
	rows = 0
	for i in tp:
		val[rows,0] = i
		val[rows,1] = mean(mat[l*(rows):l*(rows+1),2])	# calc average over all neurons for each unique time point	
		rows +=1
	return val

def mean_c(mat):					
	'''Computes conductance per unique time point averaged over all neurons.
	mat has 3 columns, 0: id, 1: time-point, 2:conductance-value 
	'''
	tp = unique(mat[:,1])			# get unique time points
	ids = unique(mat[:,0])			# get unique ids
	l = size(ids)	
	val = empty((len(tp),3))		# allocate memory for each time point
	rows = 0
	for i in tp:
		val[rows,0] = i
		val[rows,1] = mean(mat[l*(rows):l*(rows+1),2])	# calc exc-cond average over all neurons for each unique time point
		val[rows,2] = mean(mat[l*(rows):l*(rows+1),3])	# calc inh-cond average over all neurons for each unique time point		
		rows +=1
	return val

def mean_v_neuron(mat):					
	'''Computes time-averaged voltage per neuron.
	mat has 3 columns, 0: id, 1: time-point, 2:voltage-value 
	'''

	ids = unique(mat[:,0])			# get unique ids
	val = empty(size(ids))			# allocate memory for each neuron

	for i in arange(size(ids)):
		idx = ids == ids[i]
		val[i] = mean(mat[idx,2])	# calc time-averaged voltage for each neuron

	return val

def mean_rate(spikes,nr,ct,start,opt):
	'''Compute time-averaged firing rate per neuron over the whole simulation time
		spikes has 2 columns, 0: id, 1: time-points; time: period for which to calculate rate (in sec)
		nr is the total number of neurons in the population
		ct is the current simulation time (in ms !)
		start: startin gpoint from where to calculate the rate
		opt=0: take all time range into account; opt=1: exlude cs stimulation range; opt=2 only cs-range
		return is a vector with time-averaged firing rates per neuron;
		the average firing rate across all neurons and the percentage of neurons that have fired
	Caution: the mean firing rate is computed only for those neurons that have fired.'''		

	#pdb.set_trace()
	win = 80. #ms, window after CS to exclude if opt==1
	mat = spikes.copy()	# to avoid that variable spikes is changed
	mat = mat[mat[:,1]>start,:]	#choose range with starting point start
	cs_times=[]
	for i in arange(0,60.1,1):	
		cs_times = append(cs_times,100 + 2*i*100)	
	cs_times = cs_times[cs_times<=ct]	# get only spikes until current time
	time = (ct-start)/1000
	if opt==1:	#to exlude cs period
		for i in cs_times:	#set to zero all spikes that occur during cs stimulation
		    t_idx = (mat[:,1]>=i) & (mat[:,1]<=i+win)		
		    mat[t_idx,1] = 0		
		time = (ct - len(cs_times)*win)/1000		# get the net time if cs-stimulation periods are removed (in sec)

	if opt==2:	#only cs_period
		cs_mat = zeros((0,2))
		for i in cs_times:	#set to zero all spikes that occur during stimulation
		    t_idx = (mat[:,1]>=i) & (mat[:,1]<=i+win)		
		    cs_mat = concatenate((cs_mat,mat[t_idx,:]))
		time = len(cs_times)*win/1000.		# get the net time if cs-stimulation periods are removed (in sec)
	#pdb.set_trace()
	if len(mat)!=0:
		ids = unique(mat[:,0])		# get unique neuron ids
		rate = empty((len(ids),2))	# allocate memory for each time point
		for rows,id in enumerate(ids):
		    idx = mat[:,0] == id
		    rate[rows,0] = id
		    tmp = float(size(nonzero(mat[idx,1])))/time 	#get the number of spikes per neuron, divide by time to get rate
		    rate[rows,1] = round(tmp,1)	#important to keep 1 significant digit, otherwise e.g. round(0.3)=0
		idx = rate[:,1]>0
		mean_rate = mean(rate[idx,1])		# time-averaged mean rate of the neurons that have fired
		perc = float(sum(idx))/nr *100		# percentage of neurons that have fired from a total of nr neurons
	else:
		rate,mean_rate,perc = 0,0,0
	return rate,mean_rate,perc


def psd_sp(sp,nr_bins,nr_neur):
	''' Returns the power spectrum density, maximum power component and corresponding frequency 
		of the histogram of sp, where sp is a vector of spike times of a population of neurons.'''	
	h = histogram(sp,nr_bins)[0]*(1./nr_neur)
	F = fft(h-mean(h))
	F = F/len(h)		# normalize
	n = len(F)
	F = F[1:n/2+1]		# get only positive spectrum
	psd = abs(F)**2
	max_val = max(psd)
	idx = psd == max_val
	freq = find(idx)		
	return psd,max_val,freq,h

def psd(vec,res):
	''' Returns the power spectrum density, maximum power component and corresponding frequency 
		of the vector v. v could be e.g. the average membrane potential of all neurons in a population or a PSTH.
	    res = the sampling step (equals by default the simulation resolution'''
	
	#nfft = nextpow2(len(vec))
	#F = fft(vec-mean(vec),nfft)
	#F = F/len(vec)		# normalize
	#F = F[:nfft/2+1]	# get only positive spectrum
	#ff = fftfreq(nfft,res)[:nfft/2+1]
	
	F = abs(np.fft.fft(vec- np.mean(vec)))**2
	Fs = 1./(1.*0.001)
	ff2 = np.fft.fftfreq(len(vec))[0:len(vec/2)+1]
	ff = np.linspace(0,Fs/2,len(F)/2+1)
	px = ff[0:len(ff)/2+1]
	import pdb
	pdb.set_trace()
	psd = abs(F)**2
	max_val = max(psd)
	idx = psd == max_val
	max_val_freq = ff[idx]	
	return psd,ff,max_val,max_val_freq

def nextpow2(i):
    """ Find 2^n that is equal to or greater than. """ 
    n = 2
    while n < i:
	n = n * 2
    return n 


def psd_v(h,nr_neur):
	''' Returns the power spectrum density, maximum power component and corresponding frequency 
		of the vector v, which is a vector of the average membrane potential of all neurons in a population.'''
	F = fft(h-mean(h))
	F = F/len(h)		# normalize
	n = len(F)
	F = F[1:n/2+1]		# get only positive spectrum
	psd = abs(F)**2
	max_val = max(psd)
	idx = psd == max_val
	freq = find(idx)	
	return psd,max_val,freq

def pars2dict(p):			
	''' transform parameter object p to a dictionary, deletes distributions and
		replaces ParameterRange objects by arrays'''

	pars = dict(p)
	for i in dict(p):					# delete distributions if present; not to be stored in file
		if '_dist' in i:
			pars.__delitem__(i)
		if isinstance(p[i],NeuroTools.parameters.ParameterRange):	# if entry is a range (iterator) then store values as array 	
			 pars[i] = p[i]._values
	return pars
		
def dict2pars(pars):
	''' transforms dictionary pars to a parameterSet object, replaces arrays by 
		ParameterRange objects'''
		
	for i in pars:				# restore ParameterRange entries
		if isinstance(pars[i],ndarray):
			pars[i] = ParameterRange(pars[i])
	
	p = ParameterSet(pars)
	return p

def load_shelve(fname,key):
	''' Load object 'key' from file 'fname'. If key='' then a list of
		stored objects is returned.'''
	db = shelve.open(fname)
	if key=='':			# get a list with all entries in db
		mat = db.keys()	
	else:
		try:
			mat = db[key]
		except KeyError:
			print "No such key exists in database !"
	db.close()
	return mat

def save_shelve(mat,fname,key):					# save object in shelve-db
	''' Save 'mat' under name 'key' in file 'fname.'''
	db = shelve.open(fname)
	db[key] = mat
	db.close()
	
def save_pars(fname,p,key):					# save pars in a db with shelve (uses pickle)
	'''key could be the timestamp; if key=='' then 'value1' is used'''
	if key=='':
	    key = 'value1'
	pars = str(pars2dict(p))
	db = shelve.open(fname)
	db[key] = pars
	db.close()

def load_pars(fname,key):				# load pars from a db with shelve (uses pickle)

	db = shelve.open(fname)
	if key=='':			# get a list with all entries in db
		p = db.keys()	
	else:
		try:
			pars = db[key]
		except KeyError:
			print "No such key exists in database !"
			pars=[]
		pars_dict = dict2pars(eval(pars))
		p = ParameterSet(pars_dict)
	db.close()
	return p	

def create_name(index):
    '''create name from parameter index; to be used with neurotools parameterRange'''
    name=''
    for i in index:
	name = name + str(i)+'_'
    return str(name[:-1])	#ommit the last '_'

def save_hdf5(mat,fname,ar_name,gr_name):		
	''' Stores the ndarray 'mat' it the hdf5 file 'fname' under the array-name 'ar_name' in
		the group 'gr_name'. Caution: mat is transposed before stored in the table in order
		for Matlab to import it correctly'''
	
	h5file = tables.openFile(fname, mode = 'a')

	try:	
		group = h5file.root._f_getChild(gr_name)			#get group if it exists
	except tables.exceptions.NoSuchNodeError:							
		group = h5file.createGroup("/", gr_name)			#create group if it doesn't exist	
	if isinstance(mat,ndarray):
		mat = mat.T
	h5file.createArray(group, ar_name, mat)
	h5file.close()


def load_hdf5(fname,ar_name,gr_name):		
	''' Loads the ndarray 'mat' it the hdf5 file 'fname' under the array-name 'ar_name' in
		the group 'gr_name'. Caution: mat is transposed before stored in the table in order
		for Matlab to import it correctly'''
	
	h5file = tables.openFile(fname, mode = 'r')
	if gr_name != '':
		group = h5file.root._f_getChild(gr_name)			#get group if they exist
		child = group._g_loadChild(ar_name)	
		mat = child[:]
		mat = mat.T
	else:
		mat = str(h5file.root._f_listNodes())				# return a string of all goups

	h5file.close()	
	return mat


def del_gr_hdf5(fname,gr_name):
	''' Deletes a group within the hdf5 file'''
	
	h5file = tables.openFile(fname, mode = 'a')
	group = h5file.root._f_getChild(gr_name)
	group._f_remove(True)				#remove group from file
	h5file.close()

	
def ranges2list(p,label):
	''' tansforms ParameterRange objects into lists'''
	li = []
	for i in label:
		li.extend(list(p[i]._values))		
	return li


def save_res(fname,exc_sp,inh_sp,exc_v,inh_v):
	''' Save spikes and potentials to current folder'''
	
	db = shelve.open(fname)
	db['exc_sp'] = exc_sp
	db['inh_sp'] = inh_sp
	db['exc_v'] = exc_v
	db['inh_v'] = inh_v
	db.close()

def save_values(fname,names,values):
    ''' Save list of values given in 'values' with names given in 'names'''

    db = shelve.open(fname)
    pdb.set_trace()
    for i,name in enumerate(names):
	db[name] = values[i]
	print name
    db.close()	

def load_values(fname,names):
    ''' Load list of values with names given in 'names'''	
    db = shelve.open(fname)
    values=[]
    if names==[]:
	names = db.keys()
    for i,name in enumerate(names):
	values.append(db[name])
    db.close()	
    return values
	
def load_res(fname):
	''' Load spikes and potentials from current folder'''

	db = shelve.open(fname)
	exc_sp = db['exc_sp']
	inh_sp = db['inh_sp']
	exc_v = db['exc_v'] 
	inh_v = db['inh_v']
	db.close()

	return exc_sp,inh_sp,exc_v,inh_v

def rasterize(sp_times,fs,max_time,opt):
	'''Convert spike times sp_times to a binary stream. A bin is 1 when a spike has occurred within this bin and zero otherwise.
	fs is the sampling frequency in kHz, bin size = 1/fs
	CAUTION: rasterize ignores multiple spikes (clipping) that fall into the same bin. This might e.g. underestimate the firing rate
	opt=0: clipping; opt=1: return number of spikes per bin'''

	
	if sp_times.any():		# rasterize only if vector non-empty
		
		if max_time==[]:
			max_time = sp_times.max()	

		sp_times = sp_times[sp_times<=max_time]
		if len(sp_times)==0:
		  return array([])
		max_time = max_time * fs
		idx = sp_times * fs
		idx = array(idx.round(),dtype=int)

		vec = zeros((int(max_time),1))
		if opt==0:			#ignores multiple spikes
		    vec[idx-1] = 1		
		
		elif opt==1:			#returns number of spikes per bin
		    bins = unique(idx)
		    bins = concatenate((bins,array([bins[-1]+1])))
		    counts,idx2 = histogram(idx,bins)
		    idx2 = idx2[:-1]-1
		    vec[idx2] = array(counts,ndmin=2).T

		return vec.transpose()[0]
	else:					# if empty then return zeros
		return array(sp_times,dtype=int)	

def kde(spikes,width,res,pop_nr,opt,kernel,max_time):
	'''Perform a kernel density estimation of the spike train sp using a gaussian kernel
	   sp contains the spike times in msec, width the width of the kernel in msec and res (ms) the resolution of sp
	   opt=0 for density, op=1 for rate.'''

	#A gaussian width of 1 is roughly equivalent to a rectangular width of 2.5.
	#(the integral has to be one, therefore integral of  rect = 1/2.5 *(ones((25,1))) = 0.4 * (ones((25,1)) is 1
	#The sampling frequency fs (in ) is adjusted so that the resolution of the kernel is taken into account.For instance,
	#a kernel with sub-millisecond resolution of w=0.1 ms requires a fs=1./w=1/0.1=10 kHz
	#Returns density or rate and sampling frequency fs'''
	# see blog for definition of kde 
	
	#pdb.set_trace()
	if rank(spikes)>1:
	    print('First argument should be a vector containing the spike times')
	    return nan
	if ~any(spikes):	# if no spikes then return arrays of nans for subsequent processing
	    result = zeros((2,))*nan
	    xx = zeros((2,))*nan
	    kernel_func = zeros((2,))*nan
	    return result,kernel_func,xx
	    
	sp = spikes.copy()
	if len(sp)==0:return array([0]),1
	nr = len(sp)
	fs=1./res 

	sp = rasterize(sp,fs,max_time,1)

	if kernel == 'normal':
	  kernel_func = lambda t: 1./(sqrt(2*pi)*width) * exp(-t**2/(2*width**2))
	elif kernel == 'box':
	  kernel_func = lambda t: 1./width * ((t>-width/2)&(t<width/2))
	elif kernel == 'exp':	#has to be normalized (I used to compute eligibility trace)
	  kernel_func = lambda t:  exp(-t/(2*width**2))
    
	tt = arange(-width*10,width*10+1*res,res)	# better for the convolution that the length of the kernel is odd and tt is 10x width
	kernel = kernel_func(tt)		
	
	pdf = numpy.convolve(sp,kernel,mode='full')/nr	# !!! important to divide by the number of spikes		
	pdf = pdf *1e3	#multiply with 1000 because all calculations were done for width and res in ms
	
	xx = res*(arange(len(pdf))-len(kernel)/2.)	# get correct x-axis
	if opt==0:
		result = pdf		
	elif opt==1:
		result = pdf*nr/pop_nr		# return  population rate
	return result,kernel_func,xx


def setWeightsDist(prj,dist):

	src = prj._sources
	if src!=[]:			# for existing connections only
		ids = unique(src)
		ids2 = ids.tolist()
		ids2.append(max(ids2)+1)
		nr = histogram(src,ids2)
		value=[]
		for i in arange(size(ids)):
			value.append({'weights':(1000*array(dist.next(nr[0][i]),ndmin=1)).tolist()})	# !!! transform to nS (that's what Nest uses)
		nest.SetConnections(ids.tolist(),prj.plasticity_name,value)

def setDelaysDist(prj,dist):
	src = prj._sources
	if src!=[]:			# for existing connections only
		ids = unique(src)
		ids2 = ids.tolist()
		ids2.append(max(ids2)+1)
		nr = histogram(src,ids2)	# counts per source = counts per target 
		value=[]
		for i in arange(size(ids)):
			value.append({'delays':array(dist.next(nr[0][i]),ndmin=1).tolist()})
		nest.SetConnections(ids.tolist(),prj.plasticity_name,value)

def setConnectionValue(prj,param,val):
	'''Set connection value. prj: projection, param: the parameter to be set, eg. 'lr', val: the value, eg. 8e-4.'''
	src = prj._sources
	if src!=[]:			# for existing connections only
		ids = unique(src)
		ids2 = ids.tolist()
		ids2.append(max(ids2)+1)
		nr = histogram(src,ids2)	# counts per source = counts per target 
		value=[]
		for i in arange(size(ids)):
			value.append({param:repeat(val,nr[0][i]).tolist()})
		nest.SetConnections(ids.tolist(),prj.plasticity_name,value)


def findCommonTargets(spikes1,spikes2,traces,rng1,rng2,pl_name):
	'''src should be a list with the ids of the source nodes'''	
	idx1 = (spikes1[:,1] > rng1[0]) & (spikes1[:,1] < rng1[1])	# get ids of pop1
	idx2 = (spikes2[:,1] > rng2[0]) & (spikes2[:,1] < rng2[1])	# get ids of pop2
	id1 = sort(array(spikes1[idx1,0],dtype=int))			
	id2 = sort(array(spikes2[idx2,0],dtype=int))
	idx_tot = traces[0:size(unique(traces[:,0])),0]			# get the ids of all neurons in traces in the given order (not sorted)

	tgt=array([],dtype=int)
	for i in id1:					# get all targets of all sources from pop1 that spikes in rng1
		tgt=append(tgt,nest.GetConnections([i],pl_name)[0]['targets'])	

	val,bins = histogram(tgt,unique(tgt))
	tgt = unique(tgt)
	max_ids = bins[val > 0.5*max(val)]
	print size(max_ids)
	# get values only for specified time range1
	ind = (traces[:,1] > rng2[0]) & (traces[:,1] < rng2[1])
	traces = traces[ind,:]



	figure()
	v = traces[:,2]
	nr = len(unique(traces[:,0]))
	v = reshape(v,[len(v)/nr,nr])
	plot(v,color='b',marker='o')
	for id in max_ids:
		idx = idx_tot==id		
		plot(v[:,idx],color='r',marker='o')
	

#def distance(x,y):
    #''' Compute the euclidean distance between x and y. Rows are the points, columns the dimensions.'''
    #'''computes (y1-x1)^2 + (x2-x1)^2'''

    ## got the code from the web; is apparently very fast
    
    #pdb.set_trace()
    #d = zeros((x.shape[0],y.shape[0]),dtype=x.dtype)
    #for i in xrange(x.shape[1]):
        #diff2 = x[:,i,None] - y[:,i]
        #diff2 **= 2
        #d += diff2
    #return sqrt(d)

def min2(a,b):		#define new min to use with vectorize
    return min(a,b)


def distance(x,y,*args):
    ''' Compute the euclidean distance between x and y. Rows are the points, columns the dimensions.
    computes min((y1-x1),grid[0]-(y1-x1))^2 + min((y2-x2),grid[1]-(y2-x2))^2. grid[0]:width, grid[1]:height
    if grid=[0,0] then it's simply a grid, if grid=[>0,>0] it's a torus'''

 # got the code from the web; is apparently very fast
    
    min_array = numpy.vectorize(min2)	#function min2 applied to arrays
    grid = args[0]
    d = zeros((x.shape[0],y.shape[0]),dtype=x.dtype)
    for i in xrange(x.shape[1]-1):	#no z-axis used thus x.shape[1]-1
        diff2 = x[:,i,None] - y[:,i]
	diff2 = min_array(diff2,grid[i]-diff2)
        diff2 **= 2
        d += diff2
    return sqrt(d)

def adj_mat(src,tgt):
	''' Calculate the (l1 x l2) adjacency matrix, where rows correspond to sources and columns to targets.
	An entry [i,j] has a value of 1 if src[i] is connected to tgt[j], otherwise 0. '''

	src = array(src)
	tgt = array(tgt)

	src_un = unique(src)
	tgt_un = unique(tgt)

	l1 = len(src_un)
	l2 = len(tgt_un)

	mat = zeros((l1,l2))*0

	for id in src_un:
		idx = src ==id
		mat[id-src_un[0],tgt[idx]-tgt_un[0]] = 1

	return mat

def adj_w_mat(ids,prj_name):
    
    ''' compute the weighted adjacency matrix;
    is slower than adj_mat because of the loop used
    ids: the id of the nodes for which to compute the weights
    prj_name: name of projection'''
    
    l1 = len(ids)
    offset = min(ids)
    mat = zeros((l1,l1))
    #pdb.set_trace()
    for ii,sid in enumerate(ids):	
	info = nest.GetConnections([sid],prj_name)
	weights = info[0]['weights']
	tgt = array(info[0]['targets']) - offset
	#print ii,sid,len(tgt),len(weights),weights	
	if any(tgt):
	    mat[ii,tgt] = weights
	
    return mat
	
	
	
	
	
def compFanoFactor(sp):
	''' Compute the fano factor for the spike trains given in sp'''
	
	ids = unique(sp[:,0])

	counts = zeros((len(ids),))
	for i in arange(len(ids)):
		counts[i] = len(sp[sp[:,0]==i,:])

	FF = var(counts)/mean(counts)
	return FF

def restorePars(directory):
  '''restore files from directory'''
  cdir = '/home/vlachos/model/net2/'
  file = ['parameters','build_net','net2_05','plot_res']
  print file
  choice=input('select files to restore (e.g. [1,1,0,1] parameters build_net net2_05 : ')
  for nr,id in enumerate(choice):
    if id==1:
      shutil.copy(directory+file[nr]+'.py.bck',cdir+file[nr]+'.py')
      print 'Restored %s.\n'%file[nr]
  

def make_tmp(path):
    str_time1 = time.strftime('%Y%m')
    str_time2 = time.strftime('%Y%m%d_%H%M%S')
    if [i for i in os.listdir(path) if i==str_time1]==[]:	# create  directory for new month if it doesn't already exist
	    os.mkdir(path + str_time1)
    tmp_dir = path + str_time1 + '/' + str_time2 +'/'
    os.mkdir(tmp_dir)
    return tmp_dir,str_time2

def findPeaks(vec,thr):

    tmp = diff(vec)
    idx = (tmp > thr)*1
    idx2 = diff(idx)
    return vec[idx2==-1]

#def comment(tmp_dir):
    #''' comment results of simulation'''
    #fname = '/home/vlachos/model/log'
    #f = open(fname,'a')
    
    #descr = raw_input('Describe simulation: ')
    #print '\n'
    #comment = raw_input('Comment results: ')
    
    #f.write('\n\n')
    #f.write(time.strftime("%Y-%m-%d %H:%M:%S"))
    #f.write('\t'+tmp_dir+'\n')
    #f.write('\n\t\t\tDescription:\t'+descr)
    #f.write('\n\t\t\tComment: \t' +comment)
    #f.close()  

def create_db():
    conn = sqlite3.connect('/home/vlachos/model/sim_log/simulations_db')
    c = conn.cursor()

    # Create table
    c.execute('''create table simulations (description text, comments text, tags text, date text, dir text)''')

    #c.execute("""insert into simulations
    #         values ('Simulate Ping','Renewal works','p5','2006-01-05','/home')""")

    # Save (commit) the changes
    conn.commit()
    # We can also close the cursor if we are done with it
    c.close()

def comment(tmp_dir):
    ''' comment on results of simulation, records can be modified easily from firefox sqlite manager'''
    
    timepoint = tmp_dir[-16:-1]
    conn = sqlite3.connect('/home/vlachos/model/sim_log/simulations_db')
    cur = conn.cursor()
    
    description = raw_input('Describe simulation: ')
    print '\n'
    tags = raw_input('tags: ')
    comment = raw_input('Comment on results: ')
    if description=='':
	print 'Empty description. Aborted.'
    else:
	#date_time = time.ctime()
	date_time = time.ctime(time.mktime(time.strptime(timepoint,'%Y%m%d_%H%M%S')))
	print date_time
	comments2 = ''
	content = [description,comment,tags,date_time,tmp_dir,comments2]
	cur.execute('insert into simulations values (?,?,?,?,?,?)',content)
	conn.commit()
    cur.close()


def del_record(pattern):

    '''deletes all records in mysql file that match pattern	'''

    conn = sqlite3.connect('/home/vlachos/model/sim_log/simulations_db')
    cur = conn.cursor()
    pattern = '"%%%s%%"'%pattern
    cur.execute('delete from simulations where date like '+pattern)
    entries = cur.fetchall()
    conn.commit()
    cur.close()

def show_log(tag):

    # if tag given then returns only those entries witch match this tag;
    
    mytemplate = Template(filename='/home/vlachos/model/sim_log/templates/template1.html')

    conn = sqlite3.connect('/home/vlachos/model/sim_log/simulations_db')
    cur = conn.cursor()
    pattern = '"%%%s%%"'%tag
    cur.execute('select * from simulations where tags like '+pattern)
    entries = cur.fetchall()
    entries.reverse()	# show last entries first

    fh = open('/home/vlachos/model/sim_log/log.html','w')
    fh.write("<html> <head>	<h1> Simulation log </h1></head><body>")

    for entry in entries:
	print >>fh,mytemplate.render(content=entry)
    
    fh.write("</body></html>")
    fh.close()
    webbrowser.open('file:////home/vlachos/model/sim_log/log.html')
    #os.system('firefox file:////home/vlachos/model/sim_log/log.html');
    cur.close()

def get_spike_counts(sp,tot_nr):
    '''return a vector with spikes counts per neuron
    and a vector with corresponding ids, sp is the array return from
    sim.getSpikes and tot_nr is the total number of the corresponding
    population'''
    
    ids = unique(sp[:,1])
    nr = len(ids)
    counts = zeros((tot_nr,))
    for i in arange(nr):
	idx = sp[:,1] == ids[i]
	counts[ids[i]] = sum(idx) 
    return counts,ids

def psth(sp,trigger,time_range,bin_size,rate_flag):
    ''' compute the psth for all neurons given in ids. sp is the array returned from
    sim.getSpikes and trigger are the timepoints of the trigger'''

    psth_tot = []
    ids = array(unique(sp[:,1]),dtype=int)	# get ids of neurons that spiked
    bins = arange(time_range[0],time_range[1]+0.1*bin_size,bin_size) # define bins
    psth_tot = zeros((len(ids),len(bins)-1))	# allocate memory

    #import pdb
    #pdb.set_trace()

    for nr,neuron_id in enumerate(ids):
	idx = sp[:,0] == neuron_id
	timepoints = sp[idx,1] 	# get spikes per neuron

	timepoints2 = []
	for trigger_nr in trigger:
	    timepoints2.extend(\
		timepoints[(timepoints >= trigger_nr + time_range[0]) & (timepoints <= trigger_nr + time_range[1])] - trigger_nr)    
    
	psth_tot[nr,:] = histogram(timepoints2,bins)[0] / float(len(trigger))
    if rate_flag==1:
	psth_tot = psth_tot /  (bin_size/1000.)
    
    return psth_tot,bins,bin_size

def psth_conductance(conductances,trigger, time_range,res,ct):
    ''' Compute CS-triggered conductances'''

    # get rid of trigger points at the borders
    trigger = trigger[(trigger>=abs(time_range[0])) & (trigger<=ct-abs(time_range[0]))]
    exc_tot = zeros((len(trigger),diff(time_range)/res))	# allocate memory
    inh_tot = zeros((len(trigger),diff(time_range)/res))	# allocate memory
    tt = arange(0,ct-res,res)
    
    for jj,trigger_id in enumerate(trigger):
	idx = (tt>=trigger_id+time_range[0]) & (tt<=trigger_id+time_range[1]-res/2.)
	exc = conductances['exc_conductance'][idx]
	inh = conductances['inh_conductance'][idx]
	exc_tot[jj,:] = exc
	inh_tot[jj,:] = inh
    
    return mean(exc_tot,0),mean(inh_tot,0)

def psth_rate(rate,trigger, time_range,res,ct):
    ''' Compute CS-triggered PSTH'''

    # get rid of trigger points at the borders
    trigger = trigger[(trigger>=abs(time_range[0])) & (trigger<=ct-abs(time_range[0]))]
    psth = zeros((len(trigger),diff(time_range)/res))	# allocate memory
    tt = arange(0,ct-res,res)
    
    for jj,trigger_id in enumerate(trigger):
	idx = (tt>=trigger_id+time_range[0]) & (tt<=trigger_id+time_range[1]-res/2.)
	exc = rate[idx]	
	psth[jj,0:len(exc)] = exc
    
    return mean(psth,0)

def zscore(mat):

    ''' return the z-scores of the matrix computed along the rows'''

    dims = shape(mat)
    mean_mat = tile(mean(mat,1),(dims[1],1)).T
    std_mat = tile(std(mat,1),(dims[1],1)).T
    std_mat[std_mat==0] = 1	# avoid dividing by zero and producing nan (see matlab)

    zmat = (mat-mean_mat)/std_mat

    return zmat

def find_peaks(timepoints,trace,thr,time_range):
    # find the peaks of trace that exceed thr
    
    if not(time_range==[]):
	trace = trace[(timepoints>=time_range[0]) & (timepoints<=time_range[1])]
	timepoints = timepoints[(timepoints>=time_range[0]) & (timepoints<=time_range[1])]
    peaks_idx = find(diff(array(diff(trace)>0,dtype=int))==-1)	
    peaks_idx = peaks_idx + 1	# add one to get correct index
    peaks_tmp = timepoints[peaks_idx]	# get time-points of peaks
    peaks = trace[peaks_idx]
    peaks_tmp = peaks_tmp[peaks>thr]
    peaks = peaks[peaks>thr]
    
    return peaks_tmp,peaks

def find_rate_peaks(rate,timepoints,std_factor,pulses_tmp,bl_range,plot_title):
    ''' find peaks and their width  in estimated rate function
        timepoints: the timepoints vector corresponding to rate
	std_factor: determines number of standard deviations from baseline rate to set threshold
	pulse_tmp: time points of pulses
	bl_range: baseline start-end '''
    #pdb.set_trace()
    idx = (timepoints>=bl_range[0]) & (timepoints<=bl_range[1])	# time points at which to estimate baseline rate
    bl_m = mean(rate[idx])	# mean of baseline rate
    bl_std = std(rate[idx])	# std of baseline rate
    bl_var = var(rate[idx])
    print 'Baseline rate, mean = %.2f, std = %.2f \n'%(bl_m,bl_std)
    thr = bl_m + std_factor*bl_std
    idx2 = rate>thr	# get all indices for which rate is above mean + 1x std

    timepoints2 = timepoints[idx2]	
    peaks_idx = find(diff(array(diff(rate[idx2])>0,dtype=int))==-1)	# get the corresponding peaks_idx
    peaks_idx = peaks_idx +1	# add one to get correct index
    peaks_tmp = timepoints2[peaks_idx]	# get time-points of peaks
    idx = peaks_tmp>20 # ignore the first peak (below 20 ms) due to rate change
    peaks_tmp = peaks_tmp[idx]	
    peaks_idx = peaks_idx[idx]
  
    # compute duration for which rate is above mean + n x std
    tmp = timepoints.copy()*0
    tmp[idx2] = 1
    up_tp = find(diff(tmp)==1)	# find the points at which the rate goes up
    down_tp = find(diff(tmp)==-1)	# find the points at which the rate goes down

    pdb.set_trace()	
    duration = down_tp-up_tp
    
    print rate[idx2][peaks_idx],duration
    print 'Threshold: %.2f'%thr

    print 'Peaks_tmp: ', peaks_tmp
    # check false-positives and false-negatives
    # find which peaks occur within 30 ms of a given pulse
    tp = zeros((len(pulses_tmp),))	# true positives
    for ii in arange(len(pulses_tmp)):
	tp[ii] = any(abs(peaks_tmp - pulses_tmp[ii])<30) * 1
    
    tp = sum(tp)		# true positives
    fn = len(pulses_tmp) - tp 	# false negatives
    fp = len(peaks_tmp) - tp	# false positives
    recall = tp / (tp+fn)	
    precision = tp / (tp+fp)
    if (precision+recall)>0:
	f_measure = 2* (precision*recall)/(precision+recall)
    else:
	f_measure = 0

    tmp2 = mean(rate[idx2][peaks_idx])	# peaks values
    MI = ((tmp2-bl_m)/bl_m)*100		# normalize by baseline mean, in %
    
    zscore = (tmp2-bl_m)/bl_var
    difference1 = tmp2 - bl_m
    
    #pdb.set_trace()
    if plot_title!='':
	
	figure()
	plot(timepoints,rate)
	plot(peaks_tmp,rate[idx2][peaks_idx],'ro')
	plot(xlim(),[thr,thr],'k')
	plot(pulses_tmp,array(pulses_tmp)*0+thr,'go')

	title(plot_title+', MI: %.2f %% , F1: %.2f   Z: %.2f   diff: %.2f '%(MI,f_measure,zscore,difference1))
	show()
	draw()

    return rate[idx2][peaks_idx],duration,thr,MI,(tp,fp,fn,f_measure),peaks_tmp


def find_trace_peaks(trace,timepoints,std_factor,pulses_tmp,bl_range,window,plot_title):
    ''' find peaks and their width  in estimated rate function
        timepoints: the timepoints vector corresponding to rate
	std_factor: determines number of standard deviations from baseline rate to set threshold
	pulse_tmp: time points of pulses
	bl_range: baseline start-end '''
    #pdb.set_trace()
    idx = (timepoints>=bl_range[0]) & (timepoints<=bl_range[1])	# time points at which to estimate baseline rate
    bl_m = mean(trace[idx])	# mean of baseline rate
    bl_std = std(trace[idx])	# std of baseline rate

    print 'Baseline rate: mean = %.2f, std = %.2f \n'%(bl_m,bl_std)
    thr = bl_m + std_factor*bl_std
    idx2 = trace>thr	# get all indices for which trace is above threshold

    trace2 = trace[idx2]
    timepoints2 = timepoints[idx2]	
    # get changes in sign of slope from positive to negative indicating presence of a peak
    peaks_idx = find(diff(array(diff(trace2)>0,dtype=int))==-1)	
    peaks_idx = peaks_idx + 1	# add one to get correct index
    peaks_tmp = timepoints2[peaks_idx]	# get time-points of peaks
    idx = peaks_tmp>20 # ignore the first peak (below 20 ms) due to rate change
    peaks_tmp = peaks_tmp[idx]	
    peaks_idx = peaks_idx[idx]

    # get only peaks that are window points apart
    idx3 = find(diff(peaks_tmp)>window)
    peak_val = []
    idx4 = [0]
    idx4.extend(idx3.tolist())
    idx4.append(len(peaks_tmp))
    
    for ii in arange(len(idx4)-1):
	tmp4 = trace2[peaks_idx[idx4[ii]:idx4[ii+1]]]
	if len(tmp4)>0:
	    peak_val.append(max(tmp4))	
	    latency = timepoints2[trace2==max(tmp4)] - pulses_tmp
	    print 'latency: ', latency
    
    tmp5 = mean(peak_val)	     # modulation index	
    MI = (tmp5-bl_m)/abs(bl_m)*100   # normalize by peak values


    # Caution: all peaks are shown but only max peaks are used for MI
    pdb.set_trace()
    print 'difference: ',mean(peak_val)-thr
    if plot_title!='':
	
	figure()
	plot(timepoints,trace)
	plot(peaks_tmp,trace2[peaks_idx],'ro')
	plot(xlim(),[thr,thr],'k')
	plot(pulses_tmp,array(pulses_tmp)*0+thr,'go')

	title(plot_title+', MI: %.2f'%(MI))
	show()
	draw()
	
    return MI

def get_cc(spikes,raster_res,time_range,bin_size,plot_title):

    '''compute the pairwise correlation coefficient for all neurons
    raster_res: resolution of rasterized spike-train, bin_size=1/raster_res
    time_range: range for which to compute correlation coefficient
    plot_title: if not empty plot the results
    bin_size: bin size for histogram of coeffs'''

    if ~any(spikes):
	return nan,nan
	
    fs = 1./raster_res	# ms resolution for spike raster

    # get only spikes within time_range

    idx = (spikes[:,0]>=time_range[0]) & (spikes[:,0]<=time_range[1])
    spikes = spikes[idx,:]

    ids = unique(spikes[:,1])	# get ids of neurons that fired
    
    raster = zeros((len(ids),int(time_range[1]*fs)))
    for nr,id in enumerate(ids):
	idx = spikes[:,1]==id
	sp_times = spikes[idx,0]
	raster[nr,:] = rasterize(sp_times,fs,time_range[1],0)

    print shape(raster)
    #cc1 = corrcoef(raster)		# WRONG
    xc1 = dot(raster,raster.T)
    xc2 = triu(xc1,1).flatten()	# keep only upper triangular values

    bins = arange(0,20,bin_size)
    res = histogram(xc2,bins)
    print bins
    if plot_title!='':
	#figure()	
	bar(res[1][:-1]+bin_size/2.,res[0],width=0.9*bin_size)
	xlim([res[1][0],res[1][-1]])
	title(plot_title+',mean_xc: %.2f  number of pairs: %d'%(mean(xc2),len(xc2)))

    return mean(xc2),max(xc2)

def comp_pop_FF(rate,timepoints,time_range):
    '''compute FF of estimated population rate)'''
    idx = (timepoints>=time_range[0]) & (timepoints<=time_range[1])
    rate2 = rate[idx]
    FF = var(rate2)/mean(rate2)

    return FF
    
    
def create_pp(amp,std,centers,seed):
    # create a pulse-packet with AMP spikes normaly distributed around centers with std standard deviation
    # note: the pulse-packet is the same for all centers
    # set seed = [] do get always the same pattern
    if seed==[]:
	nprandom.seed(1)	# use this seed in order to replicate a specific pattern
    else:
	nprandom.seed(seed)
    pp_list=[]
    rnd = nprandom.normal(0,std,amp)
    rnd = sort(rnd)
    for nr,center in enumerate(centers):
	pp_list.extend(rnd+center)
    print 'Number of pps: ',nr+1
    return(array(pp_list)),nr+1
    
def saveSpikes(exc_sp,inh1_sp,inh2_sp,cs_exc_sp,tmp_dir,name):
	db = shelve.open(tmp_dir+name)
	db['exc_sp'] = exc_sp
	db['inh1_sp'] = inh1_sp
	db['inh2_sp'] = inh2_sp
	db['cs_exc_sp'] = cs_exc_sp
	db.close()

def loadSpikes(path,name):
	
	db = shelve.open(path+'spikes')
	if name=='':
	    return db.keys()
	spikes = db[name]
	db.close()
	return spikes
	
def get_degree(src):
  
  # compute out-degree for all neurons in src
  h_src = histogram(src,unique(src))
  idx = argsort(h_src[0])
  return h_src[0][idx],h_src[1][idx]
  
  
def change_ln(col,sc):
    ''' change lightness of color rgb color col by factor sc'''
    
    import colorsys
    
    h,l,s = colorsys.rgb_to_hls(col[0],col[1],col[2])
    rgb = colorsys.hls_to_rgb(h,l*sc,s)
    print rgb
    return array(rgb)    
 