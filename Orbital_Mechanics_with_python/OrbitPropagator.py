
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

import planetary_data as pd
import tools as t



def null_perts():
    return {
    'J2':False,
    'aero':False,
    'moon gravity':False,
    'solar_gravity':False
    }


class OrbitPropagator:
    
    def __init__(self,state0,tspan,dt,coes=False,degres=True,cb=pd.earth, perts=null_perts()):
        
        if coes:
            self.r0,self.v0,_=t.coes2rv(state0,degres=degres,mu=cb['mu'])
        else:
            self.r0=state0[:3]
            self.v0=state0[3:]

            
        self.y0=self.r0.tolist()+self.v0.tolist()
        self.tspan=tspan
        self.dt=dt
        self.cb=cb
                #number of steps
        self.n_steps=int(np.ceil(self.tspan/self.dt))
        
        #initialize arrays
        self.ts=np.zeros((self.n_steps,1))
        self.ys=np.zeros((self.n_steps,6))
        self.ts[0]=0
        self.ys[0,:]=self.y0
        self.step=1
        
        #initiate solver
        self.solver=ode(self.diffy_q)
        self.solver.set_integrator('lsoda')
        self.solver.set_initial_value(self.y0,0)

        #define perturbations dictionary
        self.perts=perts


        self.propagate_orbit()
    
    def propagate_orbit(self):
   
        #propagate the orbit
        while self.solver.successful() and self.step<self.n_steps:
            self.solver.integrate(self.solver.t+self.dt)
            self.ts[self.step]=self.solver.t
            self.ys[self.step]=self.solver.y
            self.step+=1
             
        self.rs=self.ys[:,:3]
        self.vs=self.ys[:,3:]

        
    def diffy_q(self,t,y):
    
        #unpack state
        rx,ry,rz,vx,vy,vz=y
        r=np.array([rx,ry,rz])
        v=np.array([vx,vy,vz])
    
        #norma del vector r
        norm_r=np.linalg.norm(r)
        
        #two body acceleration
        a=-r*self.cb['mu']/norm_r**3

        #J2 perturbation
        if self.perts['J2']:
            z2=r[2]**2
            r2=norm_r**2
            tx=r[0]/norm_r*(5*z2/r2-1)
            ty=r[1]/norm_r*(5*z2/r2-1)
            tz=r[2]/norm_r*(5*z2/r2-3)

            a_j2=1.5*self.cb['J2']*self.cb['mu']*self.cb['radius']**2/norm_r**4*np.array([tx,ty,tz])

            a+=a_j2
       
        
        return [vx,vy,vz,a[0],a[1],a[2]]

    def calculate_coes(self, degres=True):

        print('Calculating COEs....')

        self.coes=np.zeros((self.n_steps,6))

        for n in range (self.n_steps):
            
            self.coes[n:]=t.rv2coes(self.rs[n,:],self.vs[n,:],mu=self.cb['mu'],degres=degres)





    def plot_3d(self,show_plot=False,save_plot=False,title='Test Title'):
        
        #3D plot
        fig=plt.figure(figsize=(16,8))
        ax=fig.add_subplot(111,projection='3d')
        
        #plot trayectory and starting point
        ax.plot(self.rs[:,0],self.rs[:,1],self.rs[:,2],'w',label='Trayectory')
        ax.plot([self.rs[0,0]],[self.rs[0,1]],[self.rs[0,2]],'wo',label='Initial position')
               
        #plot earth
        _u,_v=np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
        _x=self.cb['radius']*np.cos(_u)*np.sin(_v)
        _y=self.cb['radius']*np.sin(_u)*np.sin(_v)
        _z=self.cb['radius']*np.cos(_v)
        ax.plot_surface(_x,_y,_z,cmap='Blues')
        
        l=self.cb['radius']*2.0
        x,y,z=[[0,0,0],[0,0,0],[0,0,0]]
        u,v,w=[[l,0,0],[0,l,0],[0,0,l]]
        ax.quiver(x,y,z,u,v,w,color='k')
        
        #check for custom axes limits
        max_val=np.max(np.abs(self.rs))
        
        #set lables and title
        ax.set_xlim([-max_val,max_val])
        ax.set_ylim([-max_val,max_val])
        ax.set_zlim([-max_val,max_val])
        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')
        ax.set_title(title)
        
        plt.legend()

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)
            

        
            
            
