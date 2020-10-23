import numpy as np
import datetime
import matplotlib.pyplot as plt
import math as m
from mpl_toolkits.mplot3d import Axes3D
plt.style.use('dark_background')

import planetary_data as pd

d2r=np.pi/180
r2d=180/np.pi

def norm(v):
    return np.linalg.norm(v)

def plot_n_orbits(rs,labels,cb=pd.earth, show_plot=False,save_plot=False,title='Many Orbits'):
    
    #3D plot
    fig=plt.figure(figsize=(16,8))
    ax=fig.add_subplot(111,projection='3d')
    
    #plot trayectory and starting point
    
    n=0
    for r in rs:
        ax.plot(r[:,0],r[:,1],r[:,2],label=labels[n])
        #ax.plot([r[0,0]],[r[0,1]],[r[0,2]],label='Initial position')
        n+=1
    
    #plot earth
    _u,_v=np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
    _x=cb['radius']*np.cos(_u)*np.sin(_v)
    _y=cb['radius']*np.sin(_u)*np.sin(_v)
    _z=cb['radius']*np.cos(_v)
    ax.plot_surface(_x,_y,_z,cmap='Blues')
    
    l=cb['radius']*2.0
    x,y,z=[[0,0,0],[0,0,0],[0,0,0]]
    u,v,w=[[l,0,0],[0,l,0],[0,0,l]]
    ax.quiver(x,y,z,u,v,w,color='k')
    
    #check for custom axes limits
    max_val=np.max(np.abs(rs))
        
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

#convert classical orbital elements to r and v vectors
def coes2rv(coes,degres=False,mu=pd.earth['mu']):
    if degres:
        a,e,i,ta,aop,raan,date=coes
        i*=d2r
        ta*=d2r
        aop*=d2r
        raan*=d2r
    else:
        a,e,i,ta,aop,raan,date=coes

    E=ecc_anomaly([ta,e],'tae')

    r_norm=a*(1-e**2)/(1+e*np.cos(ta))

    #calculate r and v vectors in perifocal frame
    r_perif=r_norm*np.array([m.cos(ta),m.sin(ta),0])
    v_perif=m.sqrt(mu*a)/r_norm*np.array([-m.sin(E),m.cos(E)*m.sqrt(1-e**2),0])

    #rotation matrix for perifocal to ECI
    perif2eci=np.transpose(eci2perif(raan,aop,i))

    #calculate r and v in inertial frame
    r=np.dot(perif2eci,r_perif)
    v=np.dot(perif2eci,v_perif)

    return r,v,date

#convert r and v vectors to coes(classical orbital elements) 
def rv2coes(r,v,mu=pd.earth['mu'],degres=False,print_results=False):
    #norm of position vector
    r_norm=norm(r)

    #especific angular momentum
    h=np.cross(r,v)
    h_norm=norm(h)

    #inclination
    i=m.acos(h[2]/h_norm)

    #eccentricity vector
    e=((norm(v)**2-mu/r_norm)*r-np.dot(r,v)*v)/mu

    #eccentricity escalar
    e_norm=norm(e)

    #node line
    N=np.cross([0,0,1],h)
    N_norm=norm(N)

    #RAAN
    raan=m.acos(N[0]/N_norm)
    if N[1]<0: raan=2*np.pi-raan #cuadrant check

    #argument of perigee
    aop=m.acos(np.dot(N,e)/N_norm/e_norm)
    if e[2]<0: aop=2*np.pi-aop #cuadrant check

    #true anomaly
    ta=m.acos(np.dot(e,r)/e_norm/r_norm)
    if np.dot(r,v)<0: ta=2*np.pi-ta

    #semi major axis
    a=r_norm*(1+e_norm*m.cos(ta))/(1-e_norm**2)

    if print_results:
        print('a',a)
        print('e',e)
        print('i',i*r2d)
        print('RAAN',raan*r2d)
        print('AOP',aop*r2d)
        print('TA',ta*r2d)

    #convert to degres if specify
    if degres: return[a,e_norm,i*r2d,ta*r2d,aop*r2d,raan*r2d]
    else: return[a,e_norm,i,ta,aop,raan] 

#inertial to perifocal rotation matrix
def eci2perif(raan,aop,i):
    row0=[-m.sin(raan)*m.cos(i)*m.sin(aop)+m.cos(raan)*m.cos(aop),m.cos(raan)*m.cos(i)*m.sin(aop)+m.sin(raan)*m.cos(aop),m.sin(i)*m.sin(aop)]
    row1=[-m.sin(raan)*m.cos(i)*m.cos(aop)-m.cos(raan)*m.sin(aop),m.cos(raan)*m.cos(i)*m.cos(aop)-m.sin(raan)*m.sin(aop),m.sin(i)*m.cos(aop)]
    row2=[m.sin(raan)*m.sin(i),-m.cos(raan)*m.sin(i),m.cos(i)]
    return np.array([row0,row1,row2])
    
#calculate excentricity anomaly
def ecc_anomaly(arr,method,tol=1e-8):
    if method=='newton':
        #newton's method for iteratively finding E
        Me,e=arr
        if Me<np.pi/2.0: E0=Me+e/2.0
        else: E0=Me-e
        for n in range(200): #arbitrary max numberof steps
            ratio=(E0-e*np.sin(E0)-Me)/(1-e*np.cos(E0));
            if abs(ratio)<tol:
                if n==0: return E0
                else: return E1
            else:
                E1=E0-ratio
                E0=E1
            #did not converge
            return False
    elif method=='tae':
        ta,e=arr
        return 2*m.atan(m.sqrt((1-e)/(1+e))*m.tan(ta/2.0))
    else:
        print('Invalid method for eccentric anomaly')

#takes text file containing TLEs and return classical orbital elements and other parameters  
def tle2coes(tle_filename,mu=pd.earth['mu']):
    #read tle file
    with open(tle_filename,'r') as f:
        lines=f.readlines()

    #separate into 3 lines
    line0=lines[0].strip() #name of satelite
    line1=lines[1].strip().split()
    line2=lines[2].strip().split()



    #epoch (year and day)
    epoch=line1[3]
    year,month,day,hour=calc_epoch(epoch)

    #collect coes

    #inclination
    i=float(line2[2])*d2r #rad
    #right ascention of ascending node
    raan=float(line2[2])*d2r #rad
    #eccentricity
    e_string=line2[4]
    e=float('0.' +e_string)
    #argument of perigee
    aop=float(line2[5])*d2r #rad
    #mean anomaly
    Me=float(line2[6])*d2r #rad
    #mean motion
    mean_motion=float(line2[7]) #revs/day
    #period
    T=1/mean_motion*24*3600 #seconds
    #semi mjor axis
    a=(T**2*mu/4.0/np.pi**2)**(1/3.0)

    #calculate eccentric anomaly
    E=ecc_anomaly([Me,e],'newton')

    #calculate the true anomaly
    ta=true_anomaly([E,e])

    return a,e,i,ta,aop,raan,[year,month,day,hour]

def calc_epoch(epoch):
    #year
    year=int('20'+epoch[:2])

    epoch=epoch[2:].split('.')

    #day of the year
    day_of_year=int(epoch[0])-1

    #decimal hour of day
    hour=float('0.'+epoch[1])*24.0

    #get year-month-day
    date=datetime.date(year,1,1)+datetime.timedelta(day_of_year)

    #exact month and date
    month=float(date.month)
    day=float(date.day)

    return year,month,day,hour

def true_anomaly(arr):
    E,e=arr
    return 2*np.arctan(np.sqrt((1+e)/(1-e)))*np.tan(E/2.0)

def tle2rv(tle_filename):
    return coes2rv(tle2coes(tle_filename))




        
    
        