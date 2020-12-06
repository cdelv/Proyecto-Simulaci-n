import math
import numpy as np

import tools as t
import planetary_data as pd

pi=math.pi

def lamberts_problem_universal_variables(r0,r1,deltat,tm=1,mu=pd.earth['mu'],tol=1e-6,max_steps=200,psi=0,psi_u=4*pi**2,psi_l=-4*pi):
	
	#calculate square root of mu parameter
	sqrt_mu=math.sqrt(mu)

	#calculate norm of position vectors
	r0_norm=t.norm(r0)
	r1_norm=t.norm(r1)

	#calculate gamma parameter
	gamma=np.dot(r0,r1)/r0_norm/r1_norm

	#calculate beta parameter
	beta=tm*math.sqrt(1-gamma**2)

	#calculate A parameter
	A=tm*math.sqrt(r0_norm*r1_norm*(1+gamma))

	#if A=0, solution cant be calculated
	if A==0:
		return np.array([0,0,0]),np.array([0,0,0])

	#initial values of c2 and c3
	c2=0.5
	c3=1/6.0

	#counter and solved variables
	step=0
	solved=False

	#while tolerance not met and not at max steps
	for n in range(max_steps):

		#calculate B paramenter
		B=r0_norm+r1_norm+A*(psi*c3-1)/math.sqrt(c2)

		#if A and B parameters out of range 
		if A>0 and B<0:

			#increase lower psi value
			psi_l+=pi

			#recalculate B parameter
			B*=-1

		#calculate universal variable cubed
		chi3=math.sqrt(B/c2)**3

		#calculate deltat_ variable
		deltat_=(chi3*c3+A*math.sqrt(B))/sqrt_mu

		#if difference between deltat variables is with in range
		if abs(deltat-deltat_)<tol:

			#set solved variable true
			solved=True

			#break out of the for loop
			break

		#check diferent between deltat and deltat_
		if deltat_<=deltat:

			#adjust lower psi value
			psi_l=psi

		else:

			#adjust upper psi value
			psi_u=psi

		#update psi, c2 and c3 values 
		psi=(psi_u+psi_l)/2.0
		c2=C2(psi)
		c3=C3(psi)
		#print(psi)

	#check if maximun steps reached
	if not solved:

		#algorithm did not converge on a psi value
		print('Lamberts problem universal variables dint converge')
		return np.array([0,0,0])+np.array([0,0,0])

	#calculate coeficients
	f=1-B/r0_norm
	g=A*math.sqrt(B/mu)
	gdot=1-B/r1_norm

	#calculate velocity vectors
	v0=(r1-f*r0)/g
	v1=(gdot*r1-r0)/g

	return v0,v1

def C2(psi):
	return (1-math.cos(math.sqrt(psi)))/psi

def C3(psi):
	return (math.sqrt(psi)-math.sin(math.sqrt(psi)))/(psi*math.sqrt(psi))







