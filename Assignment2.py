#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 22:58:27 2022

@author: kurup
"""

import numpy as np
import pandas as pd
from math import radians, cos, sin, fmod
import datetime
from matplotlib import pyplot as plt

c = 0 
r = 8
e1 = 1
e2 = 1
z = 0
s = 360/687
path = "01_data_mars_opposition_updated.csv"

def get_opposition(path):
  data = pd.read_csv(path)
  angles = np.array(data['ZodiacIndex']*30 + data['Degree'] + data['Minute']/60 + data['Second']/(60*60))
  days = [0.0]
  for i in range(0,11):
    day_aft = datetime.datetime(data['Year'][i+1],data['Month'][i+1],data['Day'][i+1],data['Hour'][i+1],data['Minute.1'][i+1])
    day_bef = datetime.datetime(data['Year'][i],data['Month'][i],data['Day'][i],data['Hour'][i],data['Minute.1'][i])
    no_days = day_aft-day_bef
    days.append(no_days.days+no_days.seconds/(60*60*24)+days[-1])
  return np.array([angles,days]).T
oppositions=get_opposition(path)

def get_cord(e1,e2): 
  x = e1*cos(radians(fmod(e2,360)))
  y = e1*sin(radians(fmod(e2,360)))
  return x,y

def intersect_point(ex,ey,m,cx,cy,r): ##### name of function
  tan_m = np.tan(np.radians(m))
  temp= ey - cy - ex * tan_m
  a = (1 + tan_m*tan_m) 
  b = (-2*cx) + (2*tan_m*temp)
  c = (temp*temp) + (cx*cx) - (r*r)
  x1,x2=np.roots([a,b,c])
  y1 = ey + (x1 - ex) * tan_m
  y2 = ey + (x2 - ex) * tan_m
  return [x1,x2,y1,y2]

def check_choice(angle,x1):
  temp1= (0<=angle and angle<=90) or (270<=angle and angle<=360)
  temp2= (x1<=0)
  temp3= (x1>=0)
  if (temp1 and temp3) or (not temp1 and temp2):
      return 1
  return 2

def MarsEquantModel(c,r,e1,e2,z,s,oppositions):
  points=np.array([[],[]])
  days=oppositions.T[1]
  data_found=oppositions.T[0]
  for i in range(len(data_found)):
    if(data_found[i]>180):
        data_found[i]-=360
  cCx,cCy=get_cord(1,c)
  eCx,eCy=get_cord(e1,e2+z)
  angle_calc= (days*s+z) % 360
  for k in range(12):
    temp_angle=angle_calc[k]
    x1,x2,y1,y2=intersect_point(eCx,eCy,temp_angle,cCx,cCy,r)
    temp_points=[x1,y1]
    if(check_choice(temp_angle,x1)==1):
      temp_points=[x1,y1]
    elif(check_choice(temp_angle,x1)==2):
      temp_points=[x2,y2]
    new_points=np.array([[temp_points[0]], [temp_points[1]]])
    points=np.append(points,new_points, axis=1)
  predicted=np.degrees(np.arctan2(points[1],points[0]))
  errors=predicted-data_found
  max_error=max(abs(errors))
  return errors,max_error



def bestOrbitInnerParams(r,s,oppositions):
  maxerror_ret= 1000
  for c in np.linspace(0,360,20): 
    for z in np.linspace(0,360,20):  
      for e1 in np.linspace(1,2,20): 
        for e2 in np.linspace(0,360,20):
          error_calc,maxerror_calc = MarsEquantModel(c,r,e1,e2,z,s,oppositions)
          if(maxerror_ret > maxerror_calc):
            c_ret,z_ret,e1_ret,e2_ret,errors_ret,maxerror_ret=c,z,e1,e2,error_calc,maxerror_calc
  return c_ret,e1_ret,e2_ret,z_ret,errors_ret,maxerror_ret

def bestS(r,oppositions):
  maxerror_ret = 2*360
  s_ret= 360/687
  lower_limit=(360-0.01)/687
  upper_limit=(360+0.01)/687
  divisions=5
  check_list=np.linspace(lower_limit,upper_limit,divisions)
  for s_iter in check_list:
    errors,maxError = bestOrbitInnerParams(r,s_iter,oppositions)[-2:]
    if maxerror_ret > maxError:
      maxerror_ret=maxError
      errors_ret=errors
      s_ret= s_iter
  return s_ret, errors_ret, maxerror_ret

def bestR(s,oppositions):
  maxerror_ret = 2*360 
  r_ret= 9
  lower_limit=9
  upper_limit=10
  divisions=10
  check_list=np.linspace(lower_limit,upper_limit,divisions)
  for r_iter in check_list:
    errors,maxError = bestOrbitInnerParams(r_iter,s,oppositions)[-2:]
    if maxerror_ret > maxError:
      maxerror_ret=maxError
      errors_ret=errors
      r_ret= r_iter
  return r_ret,errors_ret,maxerror_ret

def bestMarsOrbitParams(oppositions):
  s=360/687
  r,er,emax = bestR(s,oppositions) 
  s,er,emax = bestS(r,oppositions)
  
  c,e1,e2,z,errors,maxError = bestOrbitInnerParams(r,s,oppositions)
  return r,s,c,e1,e2,z,errors,maxError

r,s,c,e1,e2,z,errors,maxError=bestMarsOrbitParams(oppositions)
print(r,s,c,e1,e2,z,errors,maxError)

def sub_plot(c,e1,e2,oppositions):
    points=np.array([[],[]])
    days=oppositions.T[1]
    data_found=oppositions.T[0]
    for i in range(len(data_found)):
      if(data_found[i]>180):
          data_found[i]-=360
    cCx,cCy=get_cord(1,c)
    eCx,eCy=get_cord(e1,e2+z)
    angle_calc= (days*s+z) % 360
    for k in range(12):
      temp_angle=angle_calc[k]
      x1,x2,y1,y2=intersect_point(eCx,eCy,temp_angle,cCx,cCy,r)
      temp_points=[x1,y1]
      if(check_choice(temp_angle,x1)==1):
        temp_points=[x1,y1]
      elif(check_choice(temp_angle,x1)==2):
        temp_points=[x2,y2]
      new_points=np.array([[temp_points[0]], [temp_points[1]]])
      points=np.append(points,new_points, axis=1)
  
    return points


def plot(c,r,e1,e2,z,s,oppositions):
  plt.scatter(np.cos(np.radians(c)),np.sin(np.radians(c)),color='c')
  plt.scatter(0,0,color='r')
  
  points=sub_plot(c,e1,e2,oppositions)
  plt.scatter(points[0],points[1],color='b')
  e1,e2=get_cord(e1,e2)
  plt.scatter(e1,e2,color='k')

  
  for x,y in zip(points[0],points[1]):
      plt.arrow(0,0,2*x,2*y,linestyle='dotted',color='g')

      
      
  data_found=oppositions.T[0]
  for i in range(len(data_found)):
    if(data_found[i]>180):
        data_found[i]-=360
  plot_x=[]
  plot_y=[]
  for x in data_found:
      plot_x.append(1.5*r*np.cos(np.radians(x)))
      
  for x in data_found:    
      plot_y.append(1.5*r*np.sin(np.radians(x)))


  [plt.arrow(0,0,2*x,2*y,color='darkviolet') for x,y in zip(plot_x,plot_y)] 
  plt.ylim(-12,12)
  plt.xlim(-12,12)
  rads = np.radians(np.arange(0,360))
  plt.plot(np.cos(np.radians(c))+r*np.cos(rads),np.sin(np.radians(c))+r*np.sin(rads))
  plt.show()
