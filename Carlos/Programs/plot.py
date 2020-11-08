import numpy as np
import plotingfunctions as t
import planetary_data as pd
cb=pd.EARTH
if __name__ == '__main__':
	data0=np.loadtxt('OP.dat',delimiter='	')
	data1=np.loadtxt('OP1.dat',delimiter='	')
t.plot_n_orbits([data0,data1,],labels=['iss','HJ-2A',],cb=pd.EARTH,show_plot=True,save_plot=False,title='Many Orbits')
