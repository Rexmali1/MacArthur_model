# This program try to simulate a Mc'Arthur ecological model.
# The program have N species with a population ni whose depend of M resources ru.
# Mauricio Silva Tovar Noviember 14, 2022.

#Import necessary packagings
from numpy import *
import matplotlib.pyplot as plt
from math import *
from random import *
import csv
import pandas as pd
from datetime import date
from datetime import datetime

# The instruccion zeros((m,n)) create a matrix of n columns and m rows 
# The instruccion gauss(mu, sigma) create a number of 
#normal distribution. mu is the mean, and sigma is the standard deviation
#Matrix.transpose() transpose the matrix "Matrix"
#print('\n'.join([''.join(['{:4}'.format(item) for item in row]) for row in M])) print a matrix M 


def Assign_values_Gauss(M,mu, sigma): # This function assigns values to all elements of an matrix
                                      # with random values with mean mu and deviation standard sigma
  for i in range(len(M[0])):
    for j in range(len(M)):
      n=gauss(mu, sigma)
      while n<0:
        n=gauss(mu, sigma)    
      M[j,i]=n
  return M
 
def fi(ni,ru,mi,cu,t):        #This function return the time derivative of size population
  fit = cu.dot(ru)-mi
  fit = fit*ni
  return fit

def fu(ni,ru,au,ku,ci,t):
  fut = au*(ku-ru)-ci.dot(ni) #This function return the time derivative of resources levels
  fut = fut*ru
  return fut

#Set the number of species N and the number of resources M, an the heterogenity values
v   = 0.2       # Fraction between species and resources
N   = int(1)    # Number of species
M   = int(N/v)  # Number of species
c   = 1         # Mean of relationship between species and resources
s_c = 1         # Standar desviation of relationship between species and resources
k   = 10        # Mean of carrying capacity of resources 
s_k = 0         # Standar desviation of carrying capacity of resources
time = datetime.now().strftime('%d-%m-%Y, %H;%M;%S')
                
au  = 1         # Growing rate of resources
ku=zeros([M,1])
cu=zeros([N,M])
mi=ones([N,1])

Assign_values_Gauss(cu, c/N, s_c/sqrt(N))
ci=cu.transpose() 
Assign_values_Gauss(ku, k, s_k)


#Time simulation parameters
h    = 0.0015                      #Size of the intervalue to evaluate
tsim = 400                         #Time of simulation on days
ite  = int(tsim/h)                 #Number total of iterations 
ti   = linspace(0,tsim,num=ite+1)  #Array of times in the simulation

#Define the type of initial conditions
initial=int(input("Introduce a validate case: \n1) Random initial conditions \n2) Equal initial conditions \n"))
print("\n")

while initial!=1 and initial!=2:
  initial=int(input("Invalidate case, please introduce a validate case: \n1) Random initial conditions \n2) Equal initial conditions \n"))
  print("\n")

#Random initial conditions
if initial==1:
  from numpy import *
  ni=100*random.rand(N,1)
  ni=ni-ni%1
  ru=10*random.rand(M,1)
  
#Equal initial conditions
elif initial==2:
  ni=full([N,1],10)
  ru=full([M,1],4.5)
    
t=0

ni_s=zeros([N,ite+1])
ru_s=zeros([M,ite+1])

ni_s[:,0]=ni.transpose()
ru_s[:,0]=ru.transpose()

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

  ni_s[:,i+1]=ni.transpose()
  ru_s[:,i+1]=ru.transpose()


#ni_s=around(ni_s,decimals=0)

ni_m = mean(ni_s, axis=0)
ru_m = mean(ru_s, axis=0)

"""
#This section saves the matrix, cu, ku, ni, ru and the means of ni and ru as a cvs file
saves_cu = pd.DataFrame(cu)
saves_cu.to_csv('MacArthur Runge–Kutta cu '+'('+str(time)+').csv', index=False)
saves_ku = pd.DataFrame(ku)
saves_ku.to_csv('MacArthur Runge–Kutta ku '+'('+str(time)+').csv', index=False)

saves_ni = pd.DataFrame(ni_s).transpose()
saves_ni.to_csv('MacArthur Runge–Kutta ni '+'('+str(time)+').csv', index=False)
saves_ru = pd.DataFrame(ru_s).transpose()
saves_ru.to_csv('MacArthur Runge–Kutta ru '+'('+str(time)+').csv', index=False)


saves_ni_m = pd.DataFrame(ni_m)
saves_ni_m.to_csv('MacArthur Runge–Kutta ni_mean '+'('+str(time)+').csv', index=False)
saves_ru_m = pd.DataFrame(ru_m)
saves_ru_m.to_csv('MacArthur Runge–Kutta ru_mean '+'('+str(time)+').csv', index=False)
"""

#Plot the mean of population and resources
fig, axs = plt.subplots(1, 2)
axs[0].plot(ti, ni_m)
axs[0].set_title("Population average")
axs[0].set_xlabel('Time [days]')
axs[0].set_ylabel('Mean population')
axs[1].plot(ti, ru_m)
axs[1].set_title("Resources average")
axs[1].set_xlabel('Time [days]')
axs[1].set_ylabel('Mean resources')


fig.tight_layout()
