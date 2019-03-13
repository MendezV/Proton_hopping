#####argumentoos
###1 grid size in x direction (integer)
###2 tf duration of the driving
###3 initial state (integer)
###4 final state  (integer)
###5 stage of the engine in which the STA will be applied (0-400, each stage is of length 100)


###################################
#PACKAGES
###################################

import numpy as np
from scipy import linalg as la
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import sys

###################################

###################################





###############################################
##CONSTANTS
###############################################

#constants
#optimal600

cosa1=int(sys.argv[1])
points_x=cosa1
points_t=4*int(int(points_x**2 /10)/4.0) #para quesea multiplo de 4
DD=-20
hbar=1.0
m=1.0
Wi=10.0 #size of the well
d=5
cent=-0.0


# create grid
tf=float(sys.argv[2])
W=Wi/2.0
xf=W
x0=-W
dt=tf/points_t
dx=(xf-x0)/(points_x)
x=np.linspace(x0,xf,points_x) #debe tener un numero impar de elemntos para que sirva NInteg
t=np.linspace(0,tf,points_t)


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

def MicroDW(j,x,DD):
    ##divide in steps of 100
    frames=100
    if j<=frames:
        V=Vsqdouble(DD*(j/float(frames-1)),0,W/2,W/2,x,x0+W/4,xf-W/4)

    if (j<=2*frames) and (j> frames):
        V=Vsqdouble(DD,DD*((j-frames)/float(frames-1)),W/2,W/2,x,x0+W/4,xf-W/4)

    if (j<=3*frames) and (j> 2*frames):
        V=Vsqdouble(DD*(1-(j-2*frames)/float(frames-1)),DD,W/2,W/2,x,x0+W/4,xf-W/4)

    if (j<=4*frames) and (j> 3*frames):
        V=Vsqdouble(0,DD*(1-(j-3*frames)/float(frames-1)),W/2,W/2,x,x0+W/4,xf-W/4)
    return V


#############################################

#############################################






#############################################
##FUNCTIONS
#############################################

#########INTEGRATION########################
from scipy.interpolate import UnivariateSpline

def Ninteg(x,func,x0,xf,dx):
    spl = UnivariateSpline(x,func, k=5, s=0)
    intfunc=np.zeros(np.size(xf))
    for i in range(np.size(xf)):
        intfunc[i]=spl.integral(x0, xf[i])

    return intfunc

def Ninteg2(x,func,xi,xf,dx):
    spl = UnivariateSpline(x,func, k=5, s=0)
    intfunc=spl.integral(xi, xf)
    return intfunc


#############################################



#########DIFFERENTIATION########################

def Nderivat(func,x):
    spl = UnivariateSpline(x,func, k=5, s=0)
    derivativas=spl.derivative()
    return derivativas(x)

def Nderivat2(func,x):
    Nderivat(func,x)
    return  Nderivat( Nderivat(func,x),x)

#############################################




#########MATRIXMULTIPLICATION########################
def Band_mult(A,x):
	n=np.size(x)
	res=np.zeros(n,dtype=complex)
	res[0]=A[1][0]*x[0]+A[0][0]*x[1]
	for i in range(1,n-1):
		res[i]=A[2][i]*x[i-1]+A[1][i]*x[i]+A[0][i]*x[i+1]

	res[n-1]=A[2][n-1]*x[n-2]+A[1][n-1]*x[n-1]
	return res
#############################################




########TRIDIAGONAL SOLVER#############
def TDMASolve(a, b, c, d):
    nmax = np.size(d)#n in the numbers is rows
    # Modify the first-row coefficients
    c[0] =c[0]/b[0] #Division by zero risk.
    d[0] = d[0]/b[0]
    for i in range(1, nmax):
        ptemp = b[i] - (a[i] * c[i-1])
        c[i] /= ptemp
        d[i] = (d[i] - a[i] * d[i-1])/ptemp
    #Back Substitution
    xres = np.zeros(nmax,dtype=complex)
    xres[-1] = d[-1]
    for i in range(-2,-nmax-1,-1):
	    xres[i] = d[i] - c[i] * xres[i+1]
    return xres

#############################################

#############################################



##########################################################################################
#INITIAL CONDITION
##########################################################################################
n_ini=int(sys.argv[3])
n_fini=int(sys.argv[4])
stage=int(sys.argv[5])

Pot=MicroDW(stage,x,DD)
Ham=-0.5*(-2.0*np.diag(np.ones(points_x))+np.diag(np.ones(points_x-1),1)+np.diag(np.ones(points_x-1),-1))/(dx**2)+np.diag(Pot)   #first hamiltonian
values2=la.eigh(Ham)


psig=values2[1].T[n_ini-1] #ground state for first potental
psie=values2[1].T[n_fini-1] #excited state for first potental
Eg=values2[0][n_ini-1]   #ground state energy for first potental
Ee=values2[0][n_fini-1]   #excited state energy for first potental

psi1=psig/np.sqrt( Ninteg2(x,psig**2,x0,xf,dx) )
psi2=psie/np.sqrt( Ninteg2(x,psie**2,x0,xf,dx) )
#psi2=psigprime/np.sqrt( Ninteg2(x,psigprime**2,x0,xf,dx) )
E1=Eg
E2=Ee
#E2=Egprime
print(values2[0][0:10])

plt.plot(x,psi1)
plt.plot(x,psi2)
plt.title("initial and target states")
plt.show()
#######################################################################
#########################################################################

######SETTUP FOR TIME EVOLUTION##################################
A0p=-0.5*np.ones(points_x)/(dx**2)
A0p[0]=0
A1p=np.ones(points_x)/(dx**2)
A2p=-0.5*np.ones(points_x)/(dx**2)
A2p[points_x-1]=0
#######################################


print("...Initial Norm:",np.sum(np.conj(psi1)*(psi1))*dx)
print("...Initial eigen-Energy:",E1)



##########################################################################################
###ITERATIONS AND TIME EVOLUTION#####
##########################################################################################

rhos=[]
psipres=psi1
for i in range(points_t):
    tvar=t[i]
    Vc= MicroDW(0,x,DD)
    Hpsi=A1p+Vc
    '''
    A0matt=-1j*A0p*dt
    A1matt=1-1j*Hpsi*dt
    A2matt=-1j*A2p*dt
    mattmult=Band_mult([A2matt,A1matt,A0matt],psipres)

    A0mattinv=1j*A0p*dt
    A1mattinv=1+1j*Hpsi*dt
    A2mattinv=1j*A2p*dt
    psinew=TDMASolve(A0mattinv, A1mattinv, A2mattinv, mattmult)

    psipres=psinew
    rhos.append(np.real(psinew)**2+np.imag(psinew)**2)
    '''
    A0matt=-1*A0p*dt
    A1matt=1-1*Hpsi*dt
    A2matt=-1*A2p*dt
    mattmult=Band_mult([A2matt,A1matt,A0matt],psipres)

    A0mattinv=1*A0p*dt
    A1mattinv=1+1*Hpsi*dt
    A2mattinv=1*A2p*dt
    psinew=TDMASolve(A0mattinv, A1mattinv, A2mattinv, mattmult)

    psipres=psinew/np.sqrt( Ninteg2(x,psinew**2,x0,xf,dx) )
    rhos.append(np.real(psinew)**2+np.imag(psinew)**2)



##########################################################################################
###POSTPROCESSING#####
##########################################################################################


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

psirho1=np.array([polarR(psi2[j]) for j in range(points_x)]).T
psirho2=np.array([polarR(psinew[j]) for j in range(points_x)]).T
phiphase1=np.array([polarThe(psi2[j]) for j in range(points_x)]).T
phiphase2=np.array([polarThe(psinew[j]) for j in range(points_x)]).T
print(Ninteg(x,psirho1**2,x0,[xf],dx)[0])
print(Ninteg(x,psirho2**2,x0,[xf],dx)[0])
#psirho= np.sqrt(abs(   np.conj(np.array(psigrid))  *np.array(psigrid)   ))
#phiphase=np.real( 1j*np.log(   np.array(psigrid) / psirho   )  )
bounds=1
plt.imshow(rhos, interpolation='nearest',aspect='auto')
plt.title(r'$\rho(x,t)$')
plt.ylabel(r'$t/t_{f}$')
plt.xlabel(r'$x$')
plt.colorbar()
plt.show()

##plots

plt.plot(psirho1**2)
plt.title("Initial and final states norm")
plt.plot(psirho2**2)
plt.show()


##plots

plt.plot(np.real(psinew))
plt.plot(np.imag(psinew))
plt.title("final state real and imag")
plt.show()

##plots

plt.plot(phiphase1**2)
plt.title("Initial and final states phase")
plt.plot(phiphase2**2)
plt.show()
