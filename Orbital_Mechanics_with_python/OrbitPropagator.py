
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

import planetary_data as pd
import tools as _t

d2r=np.pi/180
r2d=180/np.pi

def null_perts():
    return {
    'J2':False,
    'aero':False,
    'Cd':0,
    'A':0,
    'mu':0,
    'moon gravity':False,
    'solar_gravity':False

    }


class OrbitPropagator:
    
    def __init__(self,state0,tspan,dt,coes=False,degres=True,cb=pd.earth, perts=null_perts(),mass0=0):
        
        if coes:
            self.r0,self.v0,_=_t.coes2rv(state0,degres=degres,mu=cb['mu'])
        else:
            self.r0=state0[:3]
            self.v0=state0[3:]
    
        self.y0=self.r0.tolist()+self.v0.tolist()
        self.tspan=tspan
        self.dt=dt
        self.cb=cb
        self.mass=mass0

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

    	print('Propagating orbit...')
   
       	#propagate the orbit
    	while self.solver.successful() and self.step<self.n_steps:

        	#integrate step
        	self.solver.integrate(self.solver.t+self.dt)

        	#extract values from solver
        	self.ts[self.step]=self.solver.t
        	self.ys[self.step]=self.solver.y
        	self.step+=1

       	#ectract array after the propagation
    	self.ts=self.ts[:self.step]    
    	self.rs=self.ys[:self.step,:3]
    	self.vs=self.ys[:self.step,3:]
    	self.alts=(np.linalg.norm(self.rs,axis=1)-self.cb['radius']).reshape((self.step,1))
        
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

        #aero drag perturbation
        if self.perts['aero']:
        	#calculate altitude and air density
        	z=norm_r-self.cb['radius']
        	rho=_t.calc_atmospheric_density(z)

        	#calculate motion of s/c with respect to rotating atmosphere
        	v_rel=v-np.cross(self.cb['atm_rot_vector'],r)

        	drag=-v_rel*0.5*rho*np.linalg.norm(v_rel)*self.perts['Cd']*self.perts['A']/self.mass

        	a+=drag
        
        return [vx,vy,vz,a[0],a[1],a[2]]

    def calculate_coes(self, degres=True):

        print('Calculating COEs....')


        self.coes=np.zeros((self.n_steps,6))

        for n in range (self.n_steps):

            self.coes[n:]=_t.rv2coes(self.rs[n,:],self.vs[n,:],mu=self.cb['mu'],degres=degres)

    def plot_coes(self,hours=False,days=False,show_plot=False,save_plot=False,title='COEs',figsize=(16,8)):
    	print('Ploting COEs ...')

    	#create figure and axses
    	fig,axs=plt.subplots(nrows=2,ncols=3,figsize=figsize)

    	#figure titles
    	fig.suptitle(title,fontsize=20)

    	#x axis
    	if hours: 
    		ts=self.ts/3600
    		xlabel='Time (hours)'
    	elif days:
    		ts=self.ts/3600/24
    		xlabel='Time (days)'
    	else:
    		ts=self.ts
    		xlabel='Time (seconds)'

    	#plot true anomaly
    	axs[0,0].plot(ts,self.coes[:,3])
    	axs[0,0].set_title('True anomaly vs.Time')
    	axs[0,0].grid(True)
    	axs[0,0].set_ylabel('Angle (degrees)')

    	#plot semi major axis
    	axs[1,0].plot(ts,self.coes[:,0])
    	axs[1,0].set_title('Semi Major Axis vs.Time')
    	axs[1,0].grid(True)
    	axs[1,0].set_ylabel('Semi Major Axis (km)')
    	axs[1,0].set_xlabel(xlabel)

    	#plot eccentricity
    	axs[0,1].plot(ts,self.coes[:,1])
    	axs[0,1].set_title('Eccentricity vs.Time')
    	axs[0,1].grid(True)

    	#plot argument of periage
    	axs[0,2].plot(ts,self.coes[:,4])
    	axs[0,2].set_title('Argument of Periapse vs.Time')
    	axs[0,2].grid(True)

    	#plot inclination
    	axs[1,1].plot(ts,self.coes[:,2])
    	axs[1,1].set_title('Inclination vs.Time')
    	axs[1,1].grid(True)
    	axs[1,1].set_ylabel('Angle (degrees)')
    	axs[1,1].set_xlabel(xlabel)

    	#plot RAAN
    	axs[1,2].plot(ts,self.coes[:,5])
    	axs[1,2].set_title('RAAN vs.Time')
    	axs[1,2].grid(True)
    	axs[1,2].set_xlabel(xlabel)

    	if show_plot:
    		plt.show()
    	if save_plot:
    		plt.savefig(title+'.png',dpi=300)

    def plot_3d(self,show_plot=False,save_plot=False,title='Test Title'):
        
        #3D plot
        fig=plt.figure(figsize=(16,8))
        ax=fig.add_subplot(111,projection='3d')
        
        #plot trayectory and starting point
        ax.plot(self.rs[:,0],self.rs[:,1],self.rs[:,2],'w',label='Trayectory',zorder=10)
        ax.plot([self.rs[0,0]],[self.rs[0,1]],[self.rs[0,2]],'wo',label='Initial position',zorder=10)
               
        #plot earth
        _u,_v=np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
        _x=self.cb['radius']*np.cos(_u)*np.sin(_v)
        _y=self.cb['radius']*np.sin(_u)*np.sin(_v)
        _z=self.cb['radius']*np.cos(_v)
        ax.plot_surface(_x,_y,_z,cmap='Blues',zorder=0)
        
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

    def calculate_apoapse_periapse(self):
    	self.apoapse=self.coes[:,0]*(1+self.coes[:,1])
    	self.periapse=self.coes[:,0]*(1-self.coes[:,1])

    def plot_apoapse_periapse(self,hours=False,days=False,show_plot=False,save_plot=False,title='Apoapse and Periapse',dpi=500):
    	#create figure
    	plt.figure(figsize=(20,10))

    	if hours:
    		ts=self.ts/3600.0
    		x_unit='Hours'
    	elif days:
    		ts=self.ts/3600/24
    		x_unit='Days'
    	else:
    		ts=self.ts
    		x_unit='Seconds'

    	#plot each
    	plt.plot(ts,self.apoapse,'b',label='Apoapse')
    	plt.plot(ts,self.periapse,'r',label='Periapse')

    	#labels
    	plt.xlabel('Time(%s)'% x_unit)
    	plt.ylabel('Altitude (km)')

    	#other parameters
    	plt.grid(True)
    	plt.title(title)
    	plt.legend()

    	if show_plot:
    		plt.show()
    	if save_plot:
    		plt.savefig(title+'.png',dpi=dpi)

    def plot_alts(self,show_plot=False,save_plot=False,hours=False,days=False,title='Radial Distance vs. Time',figsize=(16,8),dpi=500):
    	#create figure
    	if hours:
    		ts=self.ts/3600.0
    		x_unit='Hours'
    	elif days:
    		ts=self.ts/3600/24
    		x_unit='Days'
    	else:
    		ts=self.ts
    		x_unit='Seconds'

    	plt.figure(figsize=figsize)
    	plt.plot(ts,self.alts,'w')
    	plt.grid(True)
    	plt.xlabel('Time(%s)'% x_unit)
    	plt.ylabel('Altitude (km)')
    	plt.title(title)
    	if show_plot:
    		plt.show()
    	if save_plot:
    		plt.savefig(title+'.png',dpi=dpi)


            

        
            
            
