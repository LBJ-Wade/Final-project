#goal: write a program that can solve an integral symbolically and plot the results, then do it
#numerically and show that they give the same result. Also show that the galaxy is missing
#a large protion of mass. That actual curve is what think is dark matter.

# given the mass as a function of R, M(R), differentiate wrt R to get density as a function of R
#rho(R). Use that to get velocity as a function of R from freshman phys, v(R)= sqrt(G*M(R)/R), plt V(R) vs. R.
#show that this does not match the experimental data gathered. This indicates more mass we can't see
#DARK MATTER!!!!!!!!

#condiser a spiral galaxy approximated by a thin flat cylinder
import sympy as smp
from sympy import *
import matplotlib.pyplot as plt
import scipy as sci
from scipy.integrate import quad
from scipy import *
import numpy as np
from sympy.plotting import plot

#sybolic part of the code

def Symbolic():
#define constants, variables and functions
    G = 6.67191*10**-11 # m**3/kg*s**2 +_ 0.000009*10**-11
    R=Symbol('R')
    M=Function('M')

#define the M(R) function symbolically using sympy
    M = R
#get the density

#send it through the symbolic integrator
    #set the limits of integration here
    Rmin=0
    Rmax=10

    #this show the function we wish to plot
    RH = integrate(M,(R))
    print "Our Mass-Radius Equation is:"
    print(RH)
    #define the velocity function
    V = smp.sqrt(G*RH/R)
    print "Our Velocity Equation is:"
    print(V)
#now we have to prepare the antiderivative of the equation for plotting
#using the sympy plotting package.
    plot(V,(R,Rmin,Rmax))
#this returns a plot of the function M and its integral RH


if __name__=='__Symbolic__':
    Symbolic()


#time to use numpy to do evalutation numercially
#need to get this working, the integrator is not working. I want to put
#in a list and get an integrated list back. That way I can plot the function
#plus the integral
def Numerical():

    #define the params
    k=1.
    G= 6.67191*10**-11
    params =[k,G]

    #defien the density equation
    def rho(R):
        return k*R

    #set up the Radius space and the list for integration
    Rmin=0
    Rmax=100
    Rstep=1000
    R=np.linspace(Rmin,Rmax,Rstep)
    RH=[]
    #define the function desnity on R
    rho=rho(R)

    #inital conditions
    m0=1.
    M0=[m0]

    #Run the integrator over rho(R) to get M(R) as a list
    #use cumptrapz as the integrator
    M=cumtrapz(rho,M0,R,args='params')

    # take M(R) and solve for v(R)
    V=[]
    for i in range(len(R)):
        C=G*M[i]/R[i]
        if C<= 0:
            break
        else:
            return C
            v=np.sqrt(C)
            V.append(v)

    plt.figure
    plt.plot(V,R)
    plt.xlabel("Radius from the Galatic Center(x1000ly)")
    plt.ylabel("Velocity(km/s)")
    plt.title('Velocity of Mass in a Spiral Galaxy')
    plt.show()

if __name__=='__Numerical__':
    Numerical()


#put the two plots together to show simalarities and other attributes
# also allows us to
#add cumaltated error in the process and do error analysis
#do what luke said
Symbolic()
Numerical()
