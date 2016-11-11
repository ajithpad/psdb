import sys
sys.path.insert(0,'/users/padmanabhan/nest_latest/nest_scientific_linux/nest-2.2.2-build/lib/python2.7/site-packages/')
import nest
import numpy as np
import pylab as pl
import misc2
import cPickle as cp
import time
import os
import shutil
from datetime import datetime
import imp
import tempfile
import glob
import random
import pdb

class Simulation(object):
	    
    def __init__(self,prefix,new_pars=[],pars_file=[]):
        '''Creates and simulates a network in NEST'''
        
        #Temporary way to run without the selected pars file
        pars_file = []
        start_build_net = time.time()
        
        if pars_file==[]:	# import generic params_d file
	    import params_d_psdb
	    reload(params_d_psdb)
	    pars = params_d_psdb.Parameters(new_pars)		
	else: 			# import specific params_d file
	    fobj,pathname,description = imp.find_module(pars_file)
	    params_d_sp = imp.load_module(pars_file,fobj,pathname,description)
	    pars = params_d_sp.Parameters(new_pars)

	
	self.T_sim = pars.T_sim + pars.T_wup + pars.T_cdown
	#self.record_spikes = pars.record_spikes
	self.record_vm = pars.record_vm
	self.recorders = {}
	self.events = {'spikes':[],'vm':[]}	
	self.pars = pars
	self.pars.prefix = prefix	
	
	# INITIALIZE NETWORK -----------------------------------------------------------------------	        
        nest_path_tmp = tempfile.mktemp(prefix=pars.nest_path_tmp)
        os.mkdir(nest_path_tmp)
        nest.ResetKernel()
        shutil.rmtree(nest.GetStatus([0],'data_path')[0],ignore_errors=True)
        nest.SetStatus([0], {'resolution': pars.dt, 'print_time': pars.print_time,
        'overwrite_files':pars.owr_files, 'rng_seeds':[int(pars.rnd_seeds)],
        'data_path':nest_path_tmp})
        
        #print '\nBuilding network...'
        
        # CREATE SOURCES ----------------------------------------------------------------------------
        self.pg_exc = nest.Create('poisson_generator', 1)
        self.pg_inh = nest.Create('poisson_generator', 1)
        nest.SetStatus(self.pg_exc, {'rate': pars.pg_rate_exc, 'stop': pars.T_sim+pars.T_wup})
        nest.SetStatus(self.pg_inh, {'rate': pars.pg_rate_inh, 'stop': pars.T_sim+pars.T_wup})
                
        self.dc1_exc = nest.Create('dc_generator',1)
        nest.SetStatus(self.dc1_exc, pars.dc1_pars)
        
        self.dc2_exc = nest.Create('dc_generator',1)
        nest.SetStatus(self.dc2_exc, pars.dc2_pars)
        
         # CREATE POPULATIONS -----------------------------------------------------------------------
	#print 'Creating populations...\n'
	
	
	neurons_exc = []
	self.pops_exc = range(len(pars.N_exc))
	for ii,nr in enumerate(pars.N_exc):
	  self.pops_exc[ii] = nest.Create(pars.model_type, abs(nr))
	  neurons_exc.extend(self.pops_exc[ii])
	  
	#  set neuron parameters	for every population independently              
	for ntypes in range(len(pars.N_exc)):
	  nest.SetStatus(self.pops_exc[ntypes], pars.neuron_params_exc[ntypes])
	  
	
	if pars.rnd_dist:
	  nest.SetStatus(neurons_inh,'tau_m',pars.tau_m_rnd)
	  
	neurons_inh = []
	self.pops_inh = range(len(pars.N_inh))
	for ii,nr in enumerate(pars.N_inh):
	  self.pops_inh[ii] = nest.Create(pars.model_type, abs(nr))
	  neurons_inh.extend(self.pops_inh[ii])
	  
	#  set neuron parameters	for every population independently              
	for ntypes in range(len(pars.N_inh)):
	  nest.SetStatus(self.pops_inh[ntypes], pars.neuron_params_inh[ntypes])
	  
	
	if pars.rnd_dist:
	  nest.SetStatus(neurons_inh,'tau_m',pars.tau_m_rnd)
	  

	if pars.change_type:
	  self.time_lis = [pars.chg_time,pars.T_sim]
	  
	  
	self.pops = self.pops_exc + self.pops_inh
	
	self.pops_exc = [item for sublist in self.pops_exc for item in sublist]
	self.pops_inh = [item for sublist in self.pops_inh for item in sublist]
	    
	 # Make connections -------------------------------------------------------------------------
	#total_neu = [item for sublist in self.pops for item in sublist]
	self.pars.neurons_tot = len(self.pops_exc) + len(self.pops_inh) 

	self.pars.pops_exc = self.pops_exc
	self.pars.pops_inh = self.pops_inh
	
	nest.SetStatus(self.pops_exc, params = {'tau_m':pars.tau_m_exc,'C_m':pars.C_m_exc,'tau_syn':pars.tau_syn_ex})
	
	nest.DivergentConnect(self.pg_exc, self.pops_exc, weight = pars.J_ext, delay = pars.min_del)
 
	nest.DivergentConnect(self.pg_inh, self.pops_inh, weight = pars.J_ext, delay = pars.min_del)

         #STN connections
	num_stn_gpe = int(pars.epsilon_stn_gpe * len(self.pops_inh))
	nest.RandomDivergentConnect(self.pops_exc,self.pops_inh, num_stn_gpe, weight = pars.J_stn_gpe, delay = pars.delay_inter) 
	num_stn_stn = int(pars.epsilon_stn_stn* len(self.pops_exc))
	nest.RandomDivergentConnect(self.pops_exc,self.pops_exc, num_stn_stn, weight = pars.J_stn_stn, delay = pars.delay_intra) 

	#GPE connections
	num_gpe_gpe = int(pars.epsilon_gpe_gpe * len(self.pops_inh))
	nest.RandomDivergentConnect(self.pops_inh,self.pops_inh, num_gpe_gpe, weight = pars.J_gpe_gpe, delay = pars.delay_intra) 
	num_gpe_stn = int(pars.epsilon_gpe_stn* len(self.pops_exc))
	nest.RandomDivergentConnect(self.pops_inh,self.pops_exc, num_gpe_stn, weight = pars.J_gpe_stn, delay = pars.delay_inter)

	if pars.add_spikes == 1:
	    extra_spike = nest.Create('spike_generator',int(pars.num_gen))
	    def spike_times():
	      spt = np.random.uniform(pars.st_val+pars.beg_width,pars.st_val + pars.ex_width,1)
	      spike_times = (np.sort(np.random.uniform(spt,spt + pars.dur ,pars.num_spk))).tolist()     
	      return spike_times
	    

	    for ii in extra_spike:
	      spike_time = spike_times()
	      spike_time = np.around(spike_time,1)
	      nest.SetStatus([ii],params = {'spike_times':(np.array(spike_time)).flatten()})
	    
	    if self.pars.extra_spk == 'exc':
	      nest.RandomDivergentConnect(extra_spike,self.pops_inh,int(len(self.pops_inh)*pars.epsilon_stn_gpe),weight = pars.J_stn_gpe, delay = pars.delay_inter)
	      nest.RandomDivergentConnect(extra_spike,self.pops_exc,int(len(self.pops_exc)*pars.epsilon_stn_stn),weight = pars.J_stn_stn, delay = pars.delay_intra)
	    elif self.pars.extra_spk == 'inh':
	      nest.RandomDivergentConnect(extra_spike,self.pops_inh,int(len(self.pops_inh)*pars.epsilon_stn_gpe),weight = pars.J_gpe_stn, delay = pars.delay_inter)
	      nest.RandomDivergentConnect(extra_spike,self.pops_exc,int(len(self.pops_exc)*pars.epsilon_stn_stn),weight = pars.J_gpe_gpe, delay = pars.delay_intra)
 
	
	
#	pops_ch = self.pops_exc + self.pops_inh 
#	pops_ch = self.pops_exc + self.pops_inh 
#	pops_rd = random.sample(pops_ch,500)
	#CREATE RECORDERS----------------------------------------------------------------------------
	self.record_spikes = list(np.arange(self.pops_exc[0],self.pops_inh[-1]+1))
	#self.record_spikes_new = random.sample(self.record_spikes,1000) 
#	self.record_spikes = pops_rd
	#if self.record_spikes!= []:
	sd = nest.Create('spike_detector',1)
	nest.SetStatus(sd,{'to_file':True,'to_memory':True})
	nest.ConvergentConnect(self.record_spikes,sd)
	self.recorders['sd'] = sd
	
	
	if self.pars.record_vm != []:
	    vm = nest.Create('voltmeter',1)
	    #print 'Id of vm recorder: ',vm
	    nest.SetStatus(vm,{'withtime':True,'withgid':True,'to_file':True,'to_memory':False})
	    nest.DivergentConnect(vm,self.pars.record_vm)
	    nest.SetStatus(self.pars.record_vm,{'V_th':1000.}) # record free Vm
	    self.recorders['vm'] = vm
	    
	self.build_net_time = time.time()-start_build_net
	
    def sim(self):
	start_sim_time = time.time()
	if self.pars.change_type == 1:
	  for ii in self.time_lis:
	    if ii > self.pars.chg_time:
	      print 'changing type'
	      nest.SetStatus(self.pops[0][:int(self.pars.num_chg)], params = {'spb':self.pars.chg_spk_bs})
	      nest.Simulate(self.T_sim - self.pars.chg_time) 
	    else:
	      nest.Simulate(self.pars.chg_time)
	else:
	  nest.Simulate(self.T_sim)	
	self.sim_real_time = time.time() - start_sim_time
	if not os.path.exists(self.pars.data_path):
	    os.makedirs(self.pars.data_path)
	    
	self.save_spikes()
	self.save_vm()
	self.save_info()
	# important to remove any tmp folders after simulation
	tmp_dir = nest.GetStatus([0],'data_path')[0]	
	tmp_files = glob.glob(tmp_dir+'/*')
	for ff in tmp_files:	   
	    print not 'nfs' in ff
	    if not 'nfs' in ff:
		os.remove(ff)
		#print 'Removed file %s'%ff
	shutil.rmtree(tmp_dir,ignore_errors=True)
	#print 'Removed dir: ',tmp_dir
	    
    def get_spikes(self):
	  events = nest.GetStatus(self.recorders['sd'],'events')[0]
	  spikes = np.zeros((len(events['senders']),2))
	  spikes[:,0] = events['senders']
	  spikes[:,1] = events['times']
	  self.events['spikes'].append(spikes)
	    
    def get_vm(self):
	  events = np.loadtxt(nest.GetStatus(self.recorders['vm'])[0]['filenames'][0])
	  potentials = np.zeros((len(events),3))
	  potentials[:,0] = events[:,0]
	  potentials[:,1] = events[:,1]
	  potentials[:,2] = events[:,2]
	  
	  self.events['vm'].append(potentials)
	    
    def save_spikes(self,data_path=[]):
      '''merge spikes from all populations into one vector
      and save it as a numpy array'''
      
      if data_path == []:
	  data_path = self.pars.data_path
      if self.events['spikes']==[]:
	  self.get_spikes()
      spikes = self.events['spikes']
      spikes_tot = np.array([0,0],ndmin=2)
      fname = data_path + self.pars.prefix + '_spikes.npy'
      for ii,sp_pop in enumerate(spikes[0:]):
	  if len(sp_pop)>0:
	      spikes_tot = np.concatenate((spikes_tot,sp_pop))
      spikes_tot = spikes_tot[1:]	# ignore first row used for concatenation
      np.save(fname,spikes_tot)
      #print 'Saved spikes in %s'%fname
	
    def save_vm(self,data_path=[]):
	if data_path == []:
	    data_path = self.pars.data_path
	if self.events['vm']==[]:
	    vm = self.get_vm()
	vm = self.events['vm']
	vm_tot = vm[0]
	fname = data_path + self.pars.prefix + '_vm.npy'
	for ii,vm_pop in enumerate(vm[1:]):	    
	    vm_tot = np.concatenate((vm_tot,vm_pop))
	np.save(fname,vm_tot)
	#print 'Saved membrane potential in %sr'%fname
	
    def save_info(self,data_path=[]):
      if data_path == []:
	  data_path = self.pars.data_path	
      info1 = dict()
      pars0 = self.pars	
      #pars = dict()
      # convert to dict, because cPickle has problems with some objects
      pars = dict([(name,getattr(pars0,name))  for name in dir(pars0) if not name.startswith('_')])

      info1['pars'] = pars	
      tstamp = datetime.now()
      info1['tstamp'] = tstamp.strftime("%d-%m-%Y %H:%M:%S")
      info1['time_build'] = self.build_net_time
      info1['time_sim_real'] = self.sim_real_time

      pops = []
      for pop in self.pops:
	  cells = np.array(pop,dtype=int)
	  pops.append((cells[0],cells[-1]))
      info1['pops'] = pops
      #print '\n printing info:\n',info1,'\n'
      
      fname = data_path + self.pars.prefix + '.info'
      fh = open(fname,'w')
      cp.dump(info1,fh)	
      fh.close()
      #print 'Saved info-file in %s'%fname
      

	
