# This program try to simulate a Mc'Arthur ecological model.
# The program have N species with a population ni whose depend of M resources ru.
# Mauricio Silva Tovar Novie 14, 2022.

#Import necessary packagings
from numpy import *
import matplotlib.pyplot as plt
from math import *
from random import *
import csv
import pandas as pd
from datetime import date
from datetime import datetime
import os
 
if os.name == "posix":
   var = "clear"       
elif os.name == "ce" or os.name == "nt" or os.name == "dos":
   var = "cls"

# The instruccion zeros((m,n)) create a matrix of n columns and m rows 
# The instruccion gauss(mu, sigma) create a number of 
#normal distribution. mu is the mean, and sigma is the standard deviation
#Matrix.transpose() transpose the matrix "Matrix"
#print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in M])) print a matrix M 


def Assign_values_Gauss(M,mu, sigma): # This function assigns values to all elements of an matrix
                                      # with random values with mean mu and deviation standard sigma
  for i in range(len(M[0])):
    for j in range(len(M)):
      M[j,i]=gauss(mu, sigma)
  return M
 
def fi(ni,ru,mi,cu,t):      # This function return the time derivative of size population
  fit = cu.dot(ru)-1        
  fit = fit*ni              
  return fit                
                            
def fu(ni,ru,au,ku,ci,t):   # This function return the time derivative of resources levels
  fut = (ku-ru)-ci.dot(ni)
  fut = fut*ru
  return fut

v=0.3
for z in range(9):
    #Set the number of species N and the number of resources M, an the heterogenity values
    v    = v + 0.1       # Fraction between Species and resources
    N    = int(500)      # Number of species
    M    = int(N/v)      # Number of resources
    c    = 1             # Mean of relationship between species and resources
    s_c  = (0.1)**(1/2)  # Standar desviation of relationship between species and resources
    k    = 20.0          # Mean of the carrying capacity of resources
    s_k  = 0             # Standar desviation of the carrying capacity of resources
    time = datetime.now().strftime('%d-%m-%Y, %H;%M;%S')
    rep  = 20            # Number of independent runs 
    
    
    phi_s=zeros(rep)
    psi_s=zeros(rep)
    n_prom=zeros(rep)
    n_total=zeros(rep)
    n_pob=zeros([N,rep])
    
    for j in range(rep):
        au=1
        ku=zeros([M,1])
        cu=zeros([N,M])
        mi=ones([N,1])
    
        Assign_values_Gauss(cu, c/N, s_c/sqrt(N))
        ci=cu.transpose() 
        Assign_values_Gauss(ku, k, s_k)

        # Time simulation parameters
        h    = 0.001 #Size of the intervalue to evaluate
        tsim = 400 #time of simulation on days
        ite  = int(tsim/h)
        t    = 0
    
        #initial=int(input("Introduce a validate case: \n1) Random initial conditions \n2) Equal initial conditions \n"))
        #print("\n")
      
        initial=2
    
        while initial!=1 and initial!=2:
          initial=int(input("Invalidate case, please introduce a validate case: \n1) Random initial conditions \n2) Equal initial conditions \n"))
          print("\n")
        
        #Random initial conditions
        if initial==1:
          from numpy import *
          ni=100*random.rand(N,1)
          ni=ni-ni%1
          ru=10*random.rand(M,1)
          
        # Equal initial conditions
        elif initial==2:
          ni=full([N,1],10)
          ru=full([M,1],4.5)
    
        #Runge–Kutta method application 
        for i in range(ite):
    
          k1i = fi(ni,ru,mi,cu,t)
          k2i = fi(ni[:]+((h*k1i[:])/2),ru,mi,cu,t+(h/2))
          k3i = fi(ni[:]+((h*k2i[:])/2),ru,mi,cu,t+(h/2))
          k4i = fi(ni[:]+(h*k3i[:]),ru,mi,cu,t+h)
    
          k1u = fu(ni,ru,au,ku,ci,t)
          k2u = fu(ni,ru[:]+((h*k1u[:])/2),au,ku,ci,t+(h/2))
          k3u = fu(ni,ru[:]+((h*k2u[:])/2),au,ku,ci,t+(h/2))
          k4u = fu(ni,ru[:]+(h*k3u[:]),au,ku,ci,t+h)
    
          ni_v = ni + h*(k1i+2*k2i+2*k3i+k4i)/6
          ru_v = ru + h*(k1u+2*k2u+2*k3u+k4u)/6
    
          t+=h
          ni=ni_v
          ru=ru_v
       
        posArrCount1=0  # This variable count the survivor species
        
        for i in range(len(ni[:,0])):
          if (ni[i,0] >= 1/10000):
            posArrCount1 = posArrCount1 + 1
    
        posArrCount2=0  # This variable count the survivor resources
        for i in range(len(ru[:,    0])):
          if (ru[i,0] >= 1/10000):
            posArrCount2 = posArrCount2 + 1
    
        phi_s[j]   = posArrCount1/N # Fraction of survivor species 
        psi_s[j]   = posArrCount2/M # Fraction of survivor resources
        n_prom[j]  = mean(ni[:,0])
        n_total[j] = sum(ni[:,0])
        n_pob[:,j] = ni[:,0]
        os.system(var)
        print('Iteración:'+str(j))    
    
    
        #This section save the fraction of survivor species and resources and final populations
        saves_ni = pd.DataFrame(phi_s) 
        saves_ni.to_csv('MacArthur Runge–Kutta phi_s '+'m 1;a 1;c '+str(c)+';s_c^2 '+str(s_c*s_c)+';k '+str(k)+';s_k '+str(s_k)+';v '+str(v)+';N '+str(N)+'.csv', index=False)
        saves_ni =   pd.DataFrame(psi_s) 
        saves_ni.to_csv('MacArthur Runge–Kutta psi_s '+'m 1;a 1;c '+str(c)+';s_c^2 '+str(s_c*s_c)+';k '+str(k)+';s_k '+str(s_k)+';v '+str(v)+';N '+str(N)+'.csv', index=False)
        saves_ni =   pd.DataFrame(n_pob) 
        saves_ni.to_csv('MacArthur Runge–Kutta n_pob '+'m 1;a 1;c '+str(c)+';s_c^2 '+str(s_c*s_c)+';k '+str(k)+';s_k '+str(s_k)+';v '+str(v)+';N '+str(N)+'   .csv', index=False)
