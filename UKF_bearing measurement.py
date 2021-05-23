# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 17:03:39 2020

@author: Faray Majid Priyanka Das
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import math
import csv
from scipy.linalg import sqrtm
from statistics import mean


sen1=[]
sen2=[]
theta=[]
X_True=[]
Y_True=[]
Theta_True=[]
Uk=[]
Wk=[]
with open("sensor_odom.csv","r") as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=',')
  count=0
 
  for lines in csv_reader:
    if count==0:
      count=1
    else:
      sen1.append(float(lines[1]))
      sen2.append(float(lines[2]))
      theta.append(float(lines[3]))
with open("true_odometry.csv","r") as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for lines in csv_reader:
        if line_count == 0:
            line_count = 1
        else:
            X_True.append(float(lines[1]))
            Y_True.append(float(lines[2]))
            Theta_True.append(float(lines[3]))
            Uk.append(float(lines[4]))
            Wk.append(float(lines[5]))


Ts=0.033
F=np.array([[1, 0, 0],[0, 1, 0],[0, 0, 1]])
X_new=[]
Y_new=[]
Theta_new=[]

Mean=np.array([[-2],[-0.5],[0]])
Covar=np.identity(3) #Initial Estimate of the Covarriance Matrix

C1=np.array([[0, 0, 1]])
#tuning the parameter for Q and R
R=np.array([[0.1**2, 0, 0],[0, 0.1**2, 0],[0, 0, 0.1**2]])
Q=np.array([[0.1]])
rows, cols = (10, 10) 
arr = [[0]*cols]*rows 
avg=np.zeros(shape=(10,10))


#Defining Sigma Points
n=3

A1=0.4
k=0
Lamda=A1**2*(n-k) -n
Wm0=Lamda/(3+Lamda)
Wm=1/(2*(3+Lamda))
Wc=Wm
Wc0=(Lamda/3+Lamda)+(1-A1**2+2)


for o in range(2113):
  B=np.array([[Ts*math.cos(Mean[2]), 0],[Ts*math.sin(Mean[2]), 0],[0, Ts]])
  Zt=np.array([[theta[o]]])
 # print(Zt)
  Sigma=[]
  Alpha=sqrtm((n+Lamda)*Covar)
  Sigma.append(Mean)

  Delta=Mean+Alpha
  for i in range(3): 
    A=Delta[:,i]
    A.tolist()
    A=np.reshape(A,(3,1))
    Sigma.append(A)

  Delta=Mean-Alpha
  for i in range(3): 
    A=Delta[:,i]
    A.tolist()
    A=np.reshape(A,(3,1))
    Sigma.append(A)
  Beta=Sigma
  V=np.array([[Uk[o]],[Wk[o]]])
  Mean=[[0],[0],[0]]
  Covar=[[0,0,0],[0,0,0],[0,0,0]]

#Passing Sigma points through non linear funtion
#Reovering Mean and SD
  for i in range(7):
    if i==0:
      Sigma[i]=F@Sigma[i]+B@V
      Mean=Mean+(Wm0*Sigma[i])
    else:
      Sigma[i]=F@Sigma[i]+B@V
      Mean=Mean+(Wm*Sigma[i])
  for i in range(7):
    if i==0:
      D=(Sigma[i]-Beta[0])
      C=np.reshape(D,(1,3))
      Covar=Covar+((Wc0)*D@C)
    else:
      D=(Sigma[i]-Beta[0])
      C=np.reshape(D,(1,3))
      Covar=Covar+((Wc)*D@C)

  
  #C1.transpose().shape
  K=Covar@C1.transpose()@(np.linalg.inv(C1@Covar@C1.transpose()+Q))
  Mean_updt=Mean+K@(Zt-C1@Mean)
  Covar_updt=(np.identity(3)-K@C1)@Covar
  Mean=Mean_updt
  Covar=Covar_updt
  
  X_new.append(Mean[0])
  Y_new.append(Mean[1])
  Theta_new.append(Mean[2])
Time=np.linspace(start = 0, stop = 70.4, num = 2113)
plt.plot(Time,X_True, label='True Position')
plt.plot(Time,X_new, label='Estimated Position')
#plt.plot(Time,sen1, label='Sensor Position')
plt.ylabel('X Position Measurement') 
#naming the y axis 
plt.xlabel('Time') 
  
# giving a title to my graph 
plt.title('X - Position Comparrision') 
plt.legend()
plt.figure(figsize=(30, 30))  
# function to show the plot 
plt.show() 


plt.plot(Time,Y_True, label='True Position')
plt.plot(Time,Y_new, label='Estimated Position')
#plt.plot(Time,sen2, label='Sensor Position')
plt.ylabel('Y Position Measurement') 
#naming the y axis 
plt.xlabel('Time') 
  
# giving a title to my graph 
plt.title('Y - Position Comparrision') 
plt.legend()
plt.figure(figsize=(30, 30))  
# function to show the plot 
plt.show()

plt.plot(Time,Theta_True, label='True Position')
plt.plot(Time,Theta_new, label='Estimated Position')
plt.plot(Time,Theta_True, label='Sensor Position')
plt.ylabel('Theta Position Measurement') 
#naming the y axis 
plt.xlabel('Time') 
  
# giving a title to my graph 
plt.title('Theta - Position Comparrision') 
plt.legend()
plt.figure(figsize=(30, 30))  
# function to show the plot 
plt.show()
