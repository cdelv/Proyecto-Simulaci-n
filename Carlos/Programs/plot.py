import numpy as np
import plotingfunctions as t
import planetary_data as pd
cb=pd.EARTH
if __name__ == '__main__':
	data0=np.loadtxt('OP.dat',delimiter='	')
	data0=data0[:,:3]
	data1=np.loadtxt('OP1.dat',delimiter='	')
	data1=data1[:,:3]
	data2=np.loadtxt('OP2.dat',delimiter='	')
	data2=data2[:,:3]
t.plot_n_orbits([data0,data1,data2,],labels=['iss','HJ-2A','COSMOS',],cb=pd.EARTH,show_plot=True,save_plot=False,title='Many Orbits')
