
###################################
#PACKAGES
###################################

import numpy as np
from scipy import linalg as la
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import sys
from scipy import fftpack as fft

###################################

###################################





###############################################
##CONSTANTS
###############################################

#constants
#optimal 601

cosa1=601
points_x=cosa1
points_t=4*int(int(points_x**2 /25)/4.0) #para quesea multiplo de 4
DD=-20
steps=int(cosa1*5)
hbar=1.0
m=1.0
W=10.0 #size of the well
d=0.5
cent=-0.0


# create grid
tf=1000
W=W/2.0
xf=W
x0=-W
dt=tf/points_t
dx=(xf-x0)/(points_x)
x=np.linspace(x0,xf,points_x) #debe tener un numero impar de elemntos para que sirva NInteg
t=np.linspace(0,tf,points_t)
xprime=np.linspace(2*x0,2*xf,2*points_x)

#############################################

#############################################



############################################
#POTENTIALS
############################################


def window(xvec,d,cent):
    return 0.5*((np.sign(xvec+d/2.0-cent))-(np.sign(xvec-d/2.0-cent)))

def Vpike(A,xvec,d,cent):
    return A*((1-2*(xvec-cent)/d)*window(xvec,d/2,cent+d/4)+(1+2*(xvec-cent)/d)*window(xvec,d/2,cent-d/4))

def Vsq(A,xvec,d,cent):
    return A*window(xvec,d,cent)

def VJar(A1,A2,A3, xvec,cent):
    return A1*(xvec**4 - 0.5*(A2*(xvec)**2*(np.sign( (xvec-cent) )+1)+A3*(xvec)**2*(np.sign( -(xvec-cent) )+1)))

def Vharm(mw2,xvec,cent):
    return 0.5*mw2*(xvec-cent)**2

def Vsqdouble(A1,A3,d1,d3,xvec,cent1,cent2):
    return A1*window(xvec,d1,cent1)+A3*window(xvec,d3,cent2)





#############################################

#############################################






#############################################
##FUNCTIONS
#############################################

#########INTEGRATION########################
from scipy.interpolate import UnivariateSpline

def Ninteg(x,func,x0,xf,dx):
    spl = UnivariateSpline(x,func, k=3, s=0)
    derivativas=spl.derivative()
    intfunc=np.zeros(np.size(xf))
    for i in range(np.size(xf)):
        intfunc[i]=spl.integral(x0, xf[i])
 
    return intfunc

def Ninteg2(x,func,x0,xf,dx):
    spl = UnivariateSpline(x,func, k=3, s=0)
    intfunc=spl.integral(x0, xf)
    return intfunc

def eta(t):
    return ((t**3.0)/(tf**3.0))*(1 + 3.0*(1 - (t/tf)) + 6.0*(1 - (t/tf))**2.0)
#############################################


#########DIFFERENTIATION########################

def Nderivat(func,x):
    spl = UnivariateSpline(x,func, k=3, s=0)
    derivativas=spl.derivative()
    return derivativas(x)

def Nderivat2(func,x):
    Nderivat(func,x)
    return  Nderivat( Nderivat(func,x),x)

#############################################





#############################################
#INITIAL CONDITION
#############################################



psigrid=[]
norm=[]
onesies=np.diag(np.ones(points_x))
T=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)+np.diag(np.ones(points_x-1),-1))/(dx**2)+np.diag(Vsqdouble(DD,DD,2.5,2.5,x,-5+1.25,5-1.25))   #first hamiltonian

#T=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)+np.diag(np.ones(points_x-1),-1))/(dx**2)+np.diag( Vharm(1,x,0))   #first hamiltonian




T2=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)+np.diag(np.ones(points_x-1),-1))/(dx**2)+np.diag(Vsqdouble(0.05*DD,0.05*DD,2.5,2.5,x,-5+1.25,5-1.25))   #second hamiltonian


values2=la.eigh(T)
values=la.eigh(T2)

#############################################

#############################################


psi1=(1/np.sqrt(W))*(np.sin(np.pi*1*(x-W)/(2*W)))
psi1prime=(1/np.sqrt(W))*(np.sin(np.pi*1*(xprime-2*W)/(2*W)))
psi2=(1/np.sqrt(W))*(np.sin(np.pi*3*(x-W)/(2*W)))
psi2prime=(1/np.sqrt(W))*(np.sin(np.pi*3*(xprime-2*W)/(2*W)))
E1=0.5*(np.pi*1/(2*W))**2
E2=0.5*(np.pi*5/(2*W))**2
##plot initial states
'''
plt.plot(x,psi1prime)
plt.plot(x,psi2prime)
plt.show()
'''


###############################################################
filename=sys.argv[1]
V=np.loadtxt(filename,delimiter=',')
Vprime=np.concatenate((V,V), axis=1)
#Vprime=[np.concatenate((0.1*(x)**2,0.1*(x)**2), axis=None) for i in range(points_t)]
psigrid=[]
norm=[]
psigrid.append(psi1prime) #prime doulbles the size
onesies=np.diag(np.ones(points_x))
T=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)
        +np.diag(np.ones(points_x-1),-1))/(dx**2)



##plot potential
'''
bounds=1
plt.imshow(Vprime, interpolation='nearest',aspect='auto')
plt.title(r'$\rho(x,t)$')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$')
#plt.xticks(np.arange(0,(np.shape(V)[1]-2*bounds+1),(np.shape(V)[1]-2*bounds)/(2*int(xf))),np.arange(-int(xf),2*int(xf)+1))
#plt.yticks(np.arange(0,np.shape(V)[0],np.shape(V)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()
plt.show()
'''
#points_t=10000
for i in range(points_t):

	psipres=psigrid[i]
	Vc=Vprime[i][:]
	#Vc=0
	#freq =np.fft.fftshift(np.fft.fftfreq(points_x, d=dx))
	freq =np.fft.fftfreq(2*points_x, d=dx)
	#psinew=np.fft.fftshift(fft.ifft(np.exp(-(2*np.pi/4.0)*1j*dt*freq**2)*fft.fft(np.exp(-1j*dt*Vc)*psipres)))
	psinew=fft.ifft(np.exp((2*np.pi/4.0)*1j*dt*freq**2)*fft.fft(np.exp(-1j*dt*Vc)*psipres))
	#psi_new=np.concatenate( (psinew[int((psinew[points_x-1)/2):points_x], psinew[:int((psinew[points_x-1)/2)])  , axis=0)
	psigrid.append(psinew)


import math
def polarThe(z):
    a= z.real
    b= z.imag
    theta = math.atan2(b,a)
    return theta 

def polarR(z):
    a= z.real
    b= z.imag
    r = math.hypot(a,b)
    return r

psirho=np.array([[polarR(psigrid[i][j]) for i in range(points_t)] for j in range(2*points_x)]).T
phiphase=np.array([[polarThe(psigrid[i][j]) for i in range(points_t)] for j in range(2*points_x)]).T
print(Ninteg(xprime,psirho[0]**2,x0,[xf],dx)[0])
print(Ninteg(xprime,psirho[-1]**2,x0,[xf],dx)[0])
'''
psirho1=np.array([polarR(psigrid[0][j]) for j in range(2*points_x)]).T
psirho2=np.array([polarR(psigrid[-1][j]) for j in range(2*points_x)]).T
phiphase1=np.array([polarThe(psigrid[0][j]) for j in range(2*points_x)]).T
phiphase2=np.array([polarThe(psigrid[-1][j]) for j in range(2*points_x)]).T
print(Ninteg(xprime,psirho1**2,x0,[xf],dx)[0])
print(Ninteg(xprime,psirho2**2,x0,[xf],dx)[0])
'''

#psirho= np.sqrt(abs(   np.conj(np.array(psigrid))  *np.array(psigrid)   ))
#phiphase=np.real( 1j*np.log(   np.array(psigrid) / psirho   )  )



##plots
'''
plt.plot(psirho[-1]**2)
plt.plot(psi2**2)
plt.title("final state norm")
plt.show()
'''

'''
##plots

plt.plot(psirho1**2)
plt.title("Initial state norm")
plt.plot(psirho2**2)
plt.show()

#####plots

plt.plot(np.real(psigrid[-1]))
plt.plot(np.imag(psigrid[-1]))
plt.title("final state real and imag")
plt.show()

plt.plot(psigrid[-1])
plt.plot(psi2prime)
plt.title("final state and target")
plt.show()


#print((dx*np.trapz(np.conj(psi1)*psi2))*np.conj(dx*np.trapz(np.conj(psi1)*psi2)),"fidelity6")
plt.plot(phiphase1)
plt.plot(phiphase2)
plt.title("phase")
plt.show()
'''


#other fidelities
'''
	print(dx*np.abs(np.trapz(np.conj(psigrid[-1]*psi2)*(psigrid[-1]*psi2))),"fidelity1")
	print(Ninteg(x,np.abs(np.conj(psigrid[-1]*psi2)*(psigrid[-1]*psi2)),x0,[xf],dx)[0],"fidelity2")
	print((Ninteg(x,np.abs((psigrid[-1]*psi2)),x0,[xf],dx)[0])**2,"fidelity3")
	print((Ninteg(x,np.abs((psigrid[-1]*psi2)),x0,[xf],dx)[0]),"fidelity3prime")
	print((Ninteg(x,np.abs((psi2*psi2)),x0,[xf],dx)[0])**2,"fidelity4")
	print((Ninteg(x,np.abs((psi2*psi2)),x0,[xf],dx)[0]),"fidelity5")
	print((dx*np.trapz(np.conj(psi2)*psi2))*np.conj(dx*np.trapz(np.conj(psi2)*psi2)),"fidelity7")
	'''
#print((dx*np.trapz(np.conj(psigrid[-1])*psi2))*np.conj(dx*np.trapz(np.conj(psigrid[-1])*psi2)),"fidelity6")

#####plots

bounds=1
plt.imshow(psirho, interpolation='nearest',aspect='auto')
plt.title(r'$\rho(x,t)$')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$')
#plt.xticks(np.arange(0,(np.shape(V)[1]-2*bounds+1),(np.shape(V)[1]-2*bounds)/(2*int(xf))),np.arange(-int(xf),2*int(xf)+1))
#plt.yticks(np.arange(0,np.shape(V)[0],np.shape(V)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()
plt.show()

plt.imshow(phiphase, interpolation='nearest',aspect='auto')
plt.title(r'$\phi(x,t)$')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$')
#plt.xticks(np.arange(0,(np.shape(V)[1]-2*bounds+1),(np.shape(V)[1]-2*bounds)/(2*int(xf))),np.arange(-int(xf),2*int(xf)+1))
#plt.yticks(np.arange(0,np.shape(V)[0],np.shape(V)[0]/4.0),np.linspace(0,1,5))
plt.colorbar()
plt.show()

'''
np.savetxt('rho'+filename, psirho, delimiter=',')
np.savetxt('phi'+filename, phiphase, delimiter=',')
'''
