# This program try to simulate a Mc'Arthur ecological model.
# The program have N species with a population ni whose depend of M resources ru.
# Mauricio Silva Tovar August 14, 2022.

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


def Assign_values_Gauss(M,mu, sigma): #This function assigns values to all elements of an matrix
                                    #with random values with mean mu and deviation standard sigma
  for i in range(len(M[0])):
    for j in range(len(M)):
      n=gauss(mu, sigma)
      #while n<0:
        #n=gauss(mu, sigma)    
      M[j,i]=n
  return M
 
def gi(ni,ru,mi,cu,t):
  git=cu.dot(ru)-mi
  return git

def gu(ni,ru,au,ku,ci,t):
  gut=au*(ku-ru)-ci.dot(ni)
  return gut

""" 
#Set the number of species N and the number of resources M
N=int(input ("Introduce the number of species that you want to Simulate: "))
M=int(input("Introduce the number of resources implicated: "))
print("\n")

c=float(input("Please introduce the parameters c (c/N is the mean of the coficients c[i][u]): "))
s_c=float(input("Please introduce the parameters s_c((s_c)^2/N is the variance of the coficients c[i][u]): "))
print("\n")

k=float(input("Please introduce the parameters k (k is the mean of the coficients k[u]): "))
s_k=float(input("Please introduce the parameters s_k ((s_k)^2 is the variance of the coficients k[u]): "))
print("\n")
"""
#Set the number of species N and the number of resources M, an the heterogenity values
v=0.22
N=int(100)
M=int(N/v)
c=1
s_c=1
k=5
s_k=0
time = datetime.now().strftime('%d-%m-%Y, %H;%M;%S')

au=1
ku=zeros([M,1])
cu=zeros([N,M])
mi=ones([N,1])

Assign_values_Gauss(cu, c/N, s_c/sqrt(N))
ci=cu.transpose() 
Assign_values_Gauss(ku, k, s_k)

"""
Assign_values_Gauss(mi, m, s_m)
saves_mi = pd.DataFrame(ku)
saves_mi.to_csv('MacArthu Euler mi '+'('+str(time)+').csv', index=False)
"""

h=0.001 #Size of the intervalue to evaluate
tsim=100 #time of simulation on days
ite=int(tsim/h)
ti=linspace(0, tsim, num=ite+1)

initial=int(input("Introduce a validate case: \n1) Random initial conditions \n2) Equal initial conditions \n"))
print("\n")
while initial!=1 and initial!=2:
  initial=int(input("Invalidate case, please introduce a validate case: \n  1) Random initial conditions \n2) Equal initial conditions \n"))
  print("\n")

if initial==1:
  from numpy import *
  ni=100*random.rand(N,1)
  ni=ni-ni%1
  ru=10*random.rand(M,1)
elif initial==2:
  ni=full([N,1],10)
  ru=full([M,1],3)

t=0

ni_s=zeros([N,ite+1])
ru_s=zeros([M,ite+1])

ni_s[:,0]=ni.transpose()
ru_s[:,0]=ru.transpose()


for i in range(ite):
      
  ni_v=ni+h*gi(ni,ru,mi,cu,t)*ni    #Euler method application 
  ru_v=ru+h*gu(ni,ru,au,ku,ci,t)*ru #Euler method application 
  t+=h
  ni=ni_v
  #ni=around(ni,decimals=2)     #Round the population array
  #ni=ni-ni%1                   #Truncate the population array
  ru=ru_v
  ni_s[:,i+1]=ni.transpose()
  ru_s[:,i+1]=ru.transpose()


#ni_s=around(ni_s,decimals=0)

ni_m = mean(ni_s, axis=0)
ru_m = mean(ru_s, axis=0)
"""
#This section saves the matrix, cu, ku, ni, ru and the means of ni and ru as a cvs file
saves_cu = pd.DataFrame(cu)
saves_cu.to_csv('MacArthur Euler cu '+'('+str(time)+').csv', index=False)
saves_ku = pd.DataFrame(ku)
saves_ku.to_csv('MacArthur Euler ku '+'('+str(time)+').csv', index=False)

saves_ni = pd.DataFrame(ni_s).transpose()
saves_ni.to_csv('MacArthur Euler ni '+'('+str(time)+').csv', index=False)
saves_ru = pd.DataFrame(ru_s).transpose()
saves_ru.to_csv('MacArthur Euler ru '+'('+str(time)+').csv', index=False)
"""

saves_ni_m = pd.DataFrame(ni_m)
saves_ni_m.to_csv('MacArthur Euler ni_mean '+'('+str(time)+').csv', index=False)
saves_ru_m = pd.DataFrame(ru_m)
saves_ru_m.to_csv('MacArthur Euler ru_mean '+'('+str(time)+').csv', index=False)

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