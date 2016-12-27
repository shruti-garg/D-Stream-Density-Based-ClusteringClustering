# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 11:49:52 2016

@author: ic071160
"""

import math
import matplotlib.pyplot as plt
import csv
import datetime
import pandas as pd
import numpy as np

#############################################################################
############### Initializing the data structures and variables ##############
grid_list ={}                   ##{} is for dictionary, [] is for list
tgrid_list = {}
cluster = {}
input_list = []
class_name = 1

################## Functions ###################################
############### Taking the parameters as input from user #################
def param_gen(delta,cm,cl,p):
    
    ### EH: If N is too large show error and ask to decrease value of p. shape[1] gives the dimension of the input_list i.e. the number of columns
    if (cm >= N) :
        print "Parameter cm is out of range(valid range : 1<cm<N) !!!"
        exit()
    if ((cl <= 0) | (cl >= 1)) :
        print " Parameter cl is out of range(valid range : 0<cl<1) !!!"
        exit()
    d1 = float(((N - cm)/(N-cl)))
    d2 = float(cl/cm)
    d3 = float(max(d2, d1))
    d = float(math.log(d3, delta))
    dm = cm/(N*(1 - delta))
    dl = cl/(N*(1 - delta))
    gaptime = math.floor(d)  # TODO: Make sure gaptime  doesn't become zero # gaptime becomes 0 when cm is small and cl is large. #  d3 should be less than lambda
    #gaptime=1    
    assert gaptime != 0, "Modulo by zero (timediff % gaptime) !!!"
    print 'Gaptime :', gaptime
    return delta,dm,d1,gaptime,N,p
    
########### Reading the records from file to a list #################    
def get_den_grid(rec):
    grid_size= [400/50,240000/50]
    temp=rec
    key=(int((int)temp[0]/grid_size[0]),int((int)temp[1]/grid_size[1]))
    key=tuple(key)
    if key not in grid_list :
        grid_list[key]= {'tg':0,'tm':0, 'density':0,'status':'NORMAL'}                                              ##   adds key to the grid with empty values
        tgrid_list[key]= list(rec)
    else :
        tgrid_list[key].append(list(rec))
    return key
    
############### Updating the characterstic vector #######################
def update_charvec(g,d):
    global grid_list
    grid_list[g]['density']= update_den(g,d) 
    grid_list[g]['tg']= t
    return

############### Updating the density of grids ###################    
def update_den(g,den):
    global grid_list
    d = t - int(grid_list[g]['tg'])
    prev_den = grid_list[g]['density']
    grid_den = (prev_den * float(math.pow(delta,d)))+den
    return grid_den
    
###### updates all the grids by multiplying the decay factor############## 
def update_den_grids():
    global grid_list
    for grid in grid_list :
        update_charvec(grid,0)
        
################ Finding the type of grid(dense,sparse or transitional) ############    
def g_type(vector) :
    D = vector['density']
    if ( D <= dl ) :
        gtype = "SPARSE"
    elif ( D >= dm ) :
        gtype = "DENSE"
    elif (( D <= dm ) and ( D >= dl )) :
        gtype = "TRANSITIONAL"
    return gtype
        
############## Checking if the grid is an outside grid ##################### 
def isOutside(g,c):
    global cluster
    for ng in neigh_grid(g):
        if ng not in cluster[c]:
            return True
    return False
    
############## Finding the neighbour grids ######################    
def neigh_grid(g) :
    n_grid = []
    for j in range(0,len(g)):
        temp1= list(g)
        temp2 = list(g)
        temp1[j] = temp1[j]-1
        temp1 = tuple(temp1)
        temp2[j] = temp2[j]+1
        temp2 = tuple(temp2)
        if temp1 in grid_list: 
            n_grid.append(temp1)
        if temp2 in grid_list :
            n_grid.append(temp2)
    return n_grid
    
############# Move grids from one cluster to another #####################
def move(c1, c):
    global cluster
    global grid_list
    for key in cluster[c1]:
        print "Moving grid",key,"to cluster",c,"from cluster",c1 
        grid_list[key]['label'] = c
        cluster[c].append(key)
    del cluster[c1]
    print "Cluster ", c1 , "deleted"
    return
    
################ Forming the initial clusters ###########################
def initial_clustering():
    global class_name
    global grid_list
    global cluster
    
    print "Executing initial clustering..."
    update_den_grids()
    ##assign each dense grid to a unique cluster, label all other grids as NO_CLASS
    for g in grid_list :  
        t = grid_list[g]
        if(g_type(t) =="DENSE"):
            grid_list[g]['label'] = class_name            
            if class_name in cluster:
                cluster[class_name].append(g)
            else:
                cluster[class_name]= []
                cluster[class_name].append(g)
                print 'Cluster',class_name, 'created'
                class_name +=1
        else :
            grid_list[g]['label']= 'NO_CLASS'

    for c in cluster:
        for g in cluster[c]:
            if(isOutside(g,c)):
                for h in neigh_grid(g):
                    hvector= grid_list[h]
                    ch= hvector['label']
                    if (ch!= "NO_CLASS"):
                        if(len(cluster[c])> len(cluster[ch])):
                            move(ch,c)
                        else :
                            move(c,ch)
                    elif (g_type(h) == "TRANSITIONAL" ) :
                        print "Moving grid",h,"to cluster", c
                        grid_list[h]['label']= c
                        cluster[c].append(h)
                        print "Now class of", h,"is",grid_list[h]['label']
    return
    
    
    

start_time = datetime.datetime.now()
    
########### Reading the records from file to a list #################
text_file = open("test3.csv", "r")
reader = csv.reader(text_file,delimiter = ' ')
for row in reader :
    input_list.append(row)  ## TODO: Normalize the data streams

####################Input parameters from the user###################
cm=38  
cl=0.67
p=50
delta=0.998

delta,dm,dl,gaptime,N,p=param_gen(delta,cm,cl,p)
#grid_size = [400/50,240000/50] ##############input_list.apply(lambda x:x.max()/p) #@# ##TODO: put 1 instead of x.max() if the data is normalized

tc=0
t=1
for rec in input_list: #(0,input_list.shape[0]):
   # rec=input_list.iloc[i,:]   ## iloc gets the input row at location i 
    g= get_den_grid(rec)  ## g stores the key for the input record    
    d =1000              ## TODO:constant factor to be multiplied to density 
    update_charvec(g,d) ##updates the density for the new input point
    timediff= t-tc
    if timediff == gaptime :
        print 'Cluster before initial clustering', cluster
        initial_clustering()
        print 'Cluster after initial clustering', cluster
    #if (timediff % gaptime)== 0:
        #clean_grid(grid_list)
        #class_name = adjust_cluster(grid_list,class_name)
        #print 'Cluster after adjusting', cluster
    t +=1
