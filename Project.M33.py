# -*- coding: utf-8 -*-
"""
Created on Sat May 22 12:46:22 2021

@author: Sogol
"""

import numpy as np
import math



'''
K -> 0
J -> 1
H -> 2
'''

cat1 = np.loadtxt( 'Catalogue1.txt' , dtype = str )
cat2 = np.loadtxt( 'Catal2.txt' , dtype = str )


new_cat = np.zeros((len(cat2) , 5))
new_cat_index = 0

for index,number  in enumerate(cat2[: , 2]):
    if number == '0':# this is K
        #print('here')
        new_cat[new_cat_index , 0] = float(cat2[index , 0])
        new_cat[new_cat_index , 1] = cat2[index , 3]
        new_cat[new_cat_index , 2] = 0
        i = np.where(cat1[: , 0] == cat2[index , 0])
        new_cat[new_cat_index , 3] = float(cat1[i , 7])
        new_cat[new_cat_index , 4] = float(cat1[i , 3] )        
        

        
        new_cat_index += 1
        
    elif number == '1':# this is J
        new_cat[new_cat_index , 0] = float(cat2[index , 0])
        new_cat[new_cat_index , 1] = 0
        new_cat[new_cat_index , 2] = cat2[index , 3]
        i = np.where(cat1[: , 0] == cat2[index , 0])
        new_cat[new_cat_index , 3] = float(cat1[i , 7])
        new_cat[new_cat_index , 4] = float(cat1[i , 3] )  
        

        
        new_cat_index += 1
        
    elif number == '2': # this is H
        pass
    
    #new_cat[: , 0] = cat2[: , 0]

np.savetxt('data.txt', new_cat, delimiter='   ',  fmt='%.3f' )



#calculating delta index for each star
i = 0
delta = np.zeros((len(new_cat) , 2))
observation_number = np.zeros((len(delta) , 2))
K = np.zeros((len(delta) , 2))
for i in range(len(new_cat)-1):
    delta[i , 0] = new_cat[i , 0] #reporting the ID of each star
    observation_number[i , 0] = new_cat[i , 0]
    K[i , 0] = new_cat[i , 0]
    if(new_cat[i , 1] != 0): #observed in filter k
        delta[i , 1] = new_cat[i , 1] - new_cat[i , 3]
    else:
        delta[i , 1] = new_cat[i , 2] - new_cat[i , 4]
      
          
    
    
    
#calculating K variable index
i = 0

while(i < 183437):
    sigmadelta = 0
    sigmadelta2 = 0
    n = 1
    while(delta[i , 0] == delta[i+n , 0]):
        n+=1 #The number of observations for each star    
    observation_number[i, 1] = n
    for a in range(n-1):
        sigmadelta = sigmadelta + abs(delta[i+a , 1])
        sigmadelta2 = sigmadelta2 + math.pow(delta[i+a , 1] , 2)
    if(sigmadelta2 != 0):
        K[i , 1] =  ((1/n)*sigmadelta)/math.sqrt((1/n)*sigmadelta2)
    i = i + n
    
  
P = np.zeros((len(new_cat) , 2))    
x = 0   
while(x < 183437):
    if(delta[x , 0] == delta[x+1 , 0]): #observations for the same star!
        if(new_cat[x , 1] != 0 and new_cat[x+1 , 1] !=0): #both observations were in filter K
            deltaav = (delta[x ,1] + delta[x+1 , 1]) / 2  #average of 2 deltas in filter K           
            P[x , 1] = math.pow(deltaav , 2) - 1
            P[x , 0] = delta[x , 0]
            x+=2
        elif(new_cat[x , 1] == 0 and new_cat[x+1 , 1] == 0): #both observations were in filter J    
             deltaav = (delta[x ,1] + delta[x+1 , 1]) / 2  #average of 2 deltas in filter J          
             P[x , 1] = math.pow(deltaav , 2) - 1
             P[x , 0] = delta[x , 0]
             x+=2
        else:
            P[x , 1] = delta[x , 1]*delta[x+1 , 1]
            P[x , 0] = delta[x , 0]
            x+=2
    else:
        x+=1
    
    
J = np.zeros((len(delta) , 2))
j = 0
y = 0
while(y < 183435):
    if(P[y , 0] != 0):
        if(P[y , 0] == P[y+2 , 0]): #more than 1 group of observation
            while(P[y , 0] == P[y+2 , 0]):
                j = j + (np.sign(P[y , 1])*(math.sqrt(abs(P[y , 1]))))
                y+=2
            j = j + (np.sign(P[y , 1])*(math.sqrt(abs(P[y , 1]))))
            J[y , 1] = j
            J[y , 0] = P[y , 0] #transfering star id into J matrix
            j = 0
            y+=1
        else: #only 1 group of observation
            J[y , 1] = (np.sign(P[y , 1])*(math.sqrt(abs(P[y , 1]))))
            J[y , 0] = P[y , 0]
            y+=1
    else:
        y+=1
         
        
    
    
    
    



'''
array1 = np.zeros((len(cat1) , 5))
# id from cat1
array1[: , 0] = cat1[: , 0] 
# Jmag from cat1
array1[: , 1] = cat1[: , 3]
# Kmag from cat1
array1[: , 2] = cat1[: , 7]
print ( array1 )
#arrayfloat = array.astype(np.float)
'''