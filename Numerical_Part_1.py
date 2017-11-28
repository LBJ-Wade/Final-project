import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import *
from scipy import *
import scipy as sp

def Numerical():
    alpha = 1
    beta = 1
    gamma = 1
    Rho_0 = 1
    Rs = 1
    G = 1
    params = [alpha,beta,gamma,Rho_0,Rs,G]

   
    #define R  as a linear space
    R_min = 0.1
    R_max = 60
    R_step = 10000
    R = np.linspace(R_min,R_max,R_step)
 #define the density profile
    Rho =Rho_0/ ((1+(R/Rs)**2)**(3/2))
    #Rho_0/ ( (R/Rs)**gamma * (1+(R/Rs)**alpha)**((beta-gamma)/alpha) ) Somethings worng with this one

    #Rho_0/ ((1+(R/Rs)**2)**(3/2)) works well
    
    #Get Rho as a list
    #Rho = Rho(R,params)
    
    J = 2*np.pi*R #The Jacodian for the integral
    U = J*Rho #the integrand 
    #print(U)
    
    
    #Get the mass as a list
    #requieres an integral
    M =  cumtrapz(U, x=R, dx = R_step, initial = 0)
    #M = M(Rho,R,params)
    #print M
    
    V = np.sqrt((G*M)/R)
    plt.figure
    plt.plot(R,V)
    plt.plot(R,U)
    plt.plot(R,M)
    plt.show()

if __name__ == '__Numerical__':
    Numerical()
Numerical()
