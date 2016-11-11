#Program to draw the schematic for the mechanism of how chattering neurons change the oscillation 
#stability

from numpy import *
import matplotlib.pyplot as pl
#Create random integers for spread along y-axis
labelsize = 13
ticksize = 13
dis_norm = 5.
dis_norm_dis = 13.
extra_norm = 2.

num_spk  = 300
fig = pl.figure(1, figsize = (8,12))
ax = fig.add_axes([0.55,0.7,0.4,0.275])
def osc_block(start, stop, st_height, height,col,num_spk):
      x_vals = random.normal(start,stop,num_spk)
      y_vals = random.uniform(st_height,height,num_spk)
      ax.plot(x_vals,y_vals,'.',markersize = 4.,color = col)

      

for ii in arange(950.,1061.,40.):
  # Generate the positive blocks
  osc_block(ii, dis_norm, 2000., 4000.,pl.cm.Greens(200),num_spk)
  # Generate the negative blocks
  jj = ii+25.
  osc_block(jj, dis_norm,0.,2000.,pl.cm.Blues(200),num_spk)


#add additional negative spikes

osc_block(1069.,extra_norm,500.,1000.,'r',45.)

osc_block(1083.,dis_norm_dis,2000.,4000.,pl.cm.Greens(200),150.)

osc_block(1109.,dis_norm_dis,0.,2000.,pl.cm.Blues(200),150.)

osc_block(1119.,dis_norm_dis,2000.,4000.,pl.cm.Greens(200),150.)

osc_block(1133.,dis_norm,2000.,4000.,pl.cm.Greens(200),150.)

osc_block(1158.,dis_norm,0.,2000.,pl.cm.Blues(200),150.)
ax.set_xlim(950,1165)
for tl in ax.get_xticklabels():
    tl.set_fontsize(ticksize)
    
for tl in ax.get_yticklabels():
    tl.set_fontsize(ticksize)

pl.gca().xaxis.set_major_locator(pl.NullLocator())
pl.gca().yaxis.set_major_locator(pl.NullLocator())
ax.spines['top'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('Time (ms)',fontsize = labelsize)
ax.set_xlim(1010.,1150.)
#ax.set_ylabel('Neuron',fontsize = labelsize)
ax.text(1000,4000.,'B',fontsize = labelsize + 5, style = 'normal')
pl.show()