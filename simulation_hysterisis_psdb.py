import sys
#sys.path.insert(0,'/bcfgrid/data/vlachos/nest/nest_bin/lib/python2.6/site-packages')
#sys.path.insert(0,'/bcfgrid/data/padmanabhan/nest/nest_2.1-build/lib/python2.6/site-packages')
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
	    import params_d_hysterisis
	    reload(params_d_hysterisis)
	    pars = params_d_hysterisis.Parameters(new_pars)		
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
        
        print '\nBuilding network...'
        
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
	print 'Creating populations...\n'
	
	
	neurons = []
	self.pops = range(len(pars.N))
	for ii,nr in enumerate(pars.N):
	  self.pops[ii] = nest.Create(pars.model_type, abs(nr))
	  neurons.extend(self.pops[ii])
	  
	#  set neuron parameters	for every population independently              
	for ntypes in range(len(pars.N)):
	  nest.SetStatus(self.pops[ntypes], pars.neuron_params[ntypes])
	  
	
	if pars.rnd_dist:
	  nest.SetStatus(neurons,'tau_m',pars.tau_m_rnd)
	    
	 # Make connections -------------------------------------------------------------------------
	self.exc_pop = self.pops[0]

	#Create AC generator
	self.shu_gen = nest.Create('poisson_generator',params = {'rate': self.pars.shu_input,'start':pars.st_val, 'stop':pars.st_val+7*pars.shu_width})
	
	
	J_exc = pars.J_exc 
	
	J_inh = pars.J_inh 
	
	total_neu = [item for sublist in self.pops for item in sublist]
	  
	
	nest.DivergentConnect(self.pg_exc, self.exc_pop, weight = pars.J_exc, delay = pars.min_del)
	for pop in self.pops[1:]:  
           nest.DivergentConnect(self.pg_inh, pop, weight = pars.J_exc, delay = pars.min_del)

	for pop2 in self.pops:
	    self.n_exc = int(pars.epsilon * len(pop2))
	    nest.RandomDivergentConnect(self.exc_pop, pop2, self.n_exc, weight = J_exc, delay = pars.delay)
	
	
	for pop3 in self.pops:   
	    self.n_inh = int(pars.epsilon * len(pop3))
	    nest.RandomDivergentConnect(self.pops[1], pop3, self.n_inh, weight = J_inh, delay = pars.delay)
	    if len(self.pops) > 2:
	      nest.RandomDivergentConnect(self.pops[2], pop3, self.n_inh, weight = J_inh, delay = pars.delay)

	    
	#CREATE RECORDERS----------------------------------------------------------------------------
	self.record_spikes = list(np.arange(self.pops[0][0],self.pops[-1][-1]+1))
	#if self.record_spikes!= []:
	sd = nest.Create('spike_detector',1)
	nest.SetStatus(sd,{'to_file':True,'to_memory':True})
	nest.ConvergentConnect(self.record_spikes,sd)
	self.recorders['sd'] = sd
	
	
	if self.pars.record_vm != []:
	    vm = nest.Create('voltmeter',1)
	    print 'Id of vm recorder: ',vm
	    nest.SetStatus(vm,{'withtime':True,'withgid':True,'to_file':True,'to_memory':False})
	    nest.DivergentConnect(vm,self.pars.record_vm)
	    nest.SetStatus(self.pars.record_vm,{'V_th':1000.}) # record free Vm
	    self.recorders['vm'] = vm
	    
	self.build_net_time = time.time()-start_build_net
	
    def sim(self):
	start_sim_time = time.time()
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
		print 'Removed file %s'%ff
	shutil.rmtree(tmp_dir,ignore_errors=True)
	print 'Removed dir: ',tmp_dir
	    
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
      print 'Saved spikes in %s'%fname
	
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
	print 'Saved membrane potential in %sr'%fname
	
    def save_info(self,data_path=[]):
      if data_path == []:
	  data_path = self.pars.data_path	
      info1 = dict()
      pars0 = self.pars	
      #pars = dict()
      # convert to dict, because cPickle has problems with some objects
      pars = dict([(name,getattr(pars0,name))  for name in dir(pars0) if not name.startswith('_')])
      #pars['weights_exc'] = {'name':pars['weights_exc'].name,'pars':pars['weights_exc'].parameters}
      #pars['weights_inh'] = {'name':pars['weights_inh'].name,'pars':pars['weights_inh'].parameters}
      #pars['delays'] = {'name':pars['delays'].name,'pars':pars['delays'].parameters}
      #pars['conn_method1'] = pars['conn_method1'].__dict__
      #if isinstance(pars['conn_method2'],list):
	  #tmp=[]
	  #for cmlist in pars.conn_method2:
	      #tmp.append(cmlist.__dict__)
	  #pars['conn_method2'] = tmp
      #else:
	  #pars['conn_method2'] = pars['conn_method2'].__dict__	
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
      print '\n printing info:\n',info1,'\n'
      
      fname = data_path + self.pars.prefix + '.info'
      fh = open(fname,'w')
      cp.dump(info1,fh)	
      fh.close()
      print 'Saved info-file in %s'%fname
      
      
    #def post_conn(self,pre,post,post2):
	#wa = 1
	#conn_a = []
	#f_i = pre[0]
	#l_i = pre[-1]
	#p_f_i = post[0]
	#p_l_i = post[-1]
	#for jj in range(f_i,l_i+1):
	  #del wa
	  #wa = nest.FindConnections([jj])
	  #conn_num = len(nest.GetStatus(wa))
	  #conn_node = []
	  #for i in range(conn_num):
	    #if (nest.GetStatus(wa)[i]['target']>=p_f_i) and (nest.GetStatus(wa)[i]['target']<=p_l_i):
		  #f = nest.GetStatus(wa)[i]['target']
		  #conn_node.append(f)
	  #src_tgt = (jj,conn_node)
	  #conn_a.append(src_tgt)
	#return conn_a
	
