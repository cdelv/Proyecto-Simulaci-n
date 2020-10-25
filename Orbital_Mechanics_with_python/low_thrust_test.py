
import numpy as np
import planetary_data as pd
import tools as t
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts 

#time parameters
tspan=3600*1
dt=100.0

#central body
cb=pd.earth

if __name__ == '__main__':

	perts=null_perts() #next ion engine
	perts['thrust']=0.327
	perts['isp']=4300
	perts['thrust_direction']=1

	#initial mass of space craft
	mass0=10.0 #kg

	#initial state value
	state0=[cb['radius']+1000,0.01,10.0,0.0,0.0,0.0,191031]

	op=OP(state0,tspan,dt,degres=True,coes=True,mass0=mass0,perts=perts)

	op.plot_alts(show_plot=True,hours=True)
	op.plot_3d(show_plot=True)
	op.calculate_coes()
	op.plot_coes(show_plot=True,hours=True)
	op.calculate_apoapse_periapse()
	op.plot_apoapse_periapse(show_plot=True,hours=True)