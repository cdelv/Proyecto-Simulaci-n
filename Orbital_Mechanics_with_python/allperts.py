
import numpy as np
import spiceypy as spice

import planetary_data as pd
import tools as _t
import spice_tools as st
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts


#central body
cb=pd.earth

#total time
tspan=70992.6

#tspan=38000.698

#time step
dt=0.01

date0='2020-12-06-09:02:44'


mass0=419725

if __name__ == '__main__':

        #initial conditions
        iss_coes=_t.tle2coes('iss.txt',degres=True)

        #state0=[42164.0,0.001,0.0,0.0,0.0,0.0]

        #perts dictionary
        perts=null_perts()
        perts['J2']=True

        perts['aero']=True
        perts['Cd']=2.0
        perts['A']=0.004 #km²

        perts['srp']=True
        perts['a_srp']=0.004
        perts['CR']=0.4

        #perts['BSTAR']=True
        #perts['B*']=0.000020351

        #add lunar gravity perturbation
        perts['n_bodies']=[pd.moon]+[pd.sun]+[pd.jupiter]

       
        op_ISS=OP(iss_coes,tspan,dt,coes=True,degres=True,perts=perts,date0=date0,mass0=mass0,propagator='lsoda')

        op_ISS.calculate_coes(parallel=False)

        op_ISS.printcoes()
        #op_ISS.plot_coes(days=True,show_plot=False,title='Coes ISS perturbación de la luna')

        #rs_m=np.array(op_ISS.perts['n_bodies'][0]['states'][:,:3])
        #rs_m2=np.array(op_ISS.perts['n_bodies'][1]['states'][:,:3])

        #op_ISS.plot_3d(show_plot=True)

       # _t.plot_n_orbits([op_ISS.rs],labels=['ISS'],show_plot=True)
