# D-Stream-Density-Based-Clustering
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 16:38:17 2016

@author: Z003PF5Z
"""

import math
import matplotlib.pyplot as plt
import csv
import datetime
import pandas as pd
import numpy as np

################## Functions ###################################
############### Taking the parameters as input from user #################
def param_gen(delta,cm,cl,p):
    N=p**input_list.shape[1] ### EH: If N is too large show error and ask to decrease value of p
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
    #gaptime = math.floor(d)  # TODO: Make sure gaptime  doesn't become zero # gaptime becomes 0 when cm is small and cl is large. #  d3 should be less than lambda
    gaptime=1    
    assert gaptime != 0, "Modulo by zero (timediff % gaptime) !!!"
    print 'Gaptime :', gaptime
    return delta,dm,d1,gaptime,N,p
#####################################assign or add grids for a data point####################################################
def get_den_grid(rec):
    temp=rec
    key=(temp/grid_size).apply(int)
    key=tuple(key)
    if key not in grid_list :
        grid_list[key]= {'tg':0,'tm':0, 'density':0,'status':'NORMAL'}                                              ##   adds key to the grid with empty values
        tgrid_list[key]= list(rec)
    else :
        tgrid_list[key].append(list(rec))
    return key
############### Updating the characterstic vector #######################
def update_charvec(g,d):
    vector = grid_list[g]
    vector['density']= update_den(g,d) 
    vector['tg']= t
    grid_list[g] = vector        ##i!
    return
    
############### Updating the density of grids ###################    
def update_den(g,den):
    prev_vect= grid_list[g]
    d = t - int(prev_vect['tg'])
    prev_den = prev_vect['density']
    grid_den = prev_den * float(math.pow(delta,d))+den
    return grid_den
################ Forming the initial clusters ###########################
def initial_clustering(grid_list,class_name):
    print "In Initial Clustering ...."
    print 'Clusters before', cluster
    update_den_grids(grid_list)
    print 'Grid_list while initialisation',grid_list
    for grid in grid_list :
        t = grid_list[grid]
        if(g_type(t) =="DENSE"):
            t['label'] = class_name            
            if class_name in cluster:
                cluster[class_name].append(grid)
            else:
                cluster[class_name]= []
                cluster[class_name].append(grid)
                print 'Cluster',class_name, 'created'
                class_name +=1
        else :
            print "No cluster"
            t['label']= 'NO_CLASS'
        grid_list[grid] = t   ##i!
    cluster1= dict(cluster)
    for c in cluster1:
        for g in cluster1[c]:
            if(isoutside(g,cluster1[c])):
                for hkey in neigh_grid(g):
                    h= grid_list[hkey]
                    c1= h['label']
                    if (c1!= "NO_CLASS"):
                        if(len(cluster1[c])> len(cluster1[c1])):
                            move(c1,c)
                        else :
                            move(c,c1)
                    elif (g_type(h) == "TRANSITIONAL" ) :
                        print "Moving grid",hkey,"to cluster", c
                        h['label']= c
                        cluster[c].append(hkey)
                        print "Now class of", hkey,"is",h['label']
                    grid_list[hkey] = h 
    #plot_cluster()
    return class_name
    
###### updates all the grids by multiplying the decay factor############## 
def update_den_grids(grid_list):
    for grid in grid_list :
        update_charvec(grid,0)
############## Checking if the grid is an outside grid ##################### 
def isoutside(g,ckey):
    for ng in neigh_grid(g):
        if ng not in ckey:
            return True
############# Moving a grid from one cluster to another #####################
def move(c1, c):
    for key in cluster[c1]:
        print "Moving grid",key,"to cluster",c,"from cluster",c1 
        value = grid_list[key]
        value['label']= c
        grid_list[key]=value    ##i!
        cluster[c].append(key)
    temp = c1
    del cluster[c1]
    cluster[temp] = []
################ Finding the type of grid(dense,sparse or transitional) ############    
def g_type(value) :
    D = value['density']
    if ( D <= dl ) :
        gtype = "SPARSE"
    elif ( D >= dm ) :
        gtype = "DENSE"
    elif (( D <= dm ) and ( D >= dl )) :
        gtype = "TRANSITIONAL"
    return gtype
############## Finding the neighbour grids ######################    
def neigh_grid(grid) :
    n_grid = []
    for j in range(len(grid)):
        temp1= list(grid)
        temp2 = list(grid)
        temp1[j] = temp1[j]-1
        temp1 = tuple(temp1)
        temp2[j] = temp2[j]+1
        temp2 = tuple(temp2)
        if temp1 in grid_list: 
            n_grid.append(temp1)
        if temp2 in grid_list :
            n_grid.append(temp2)
    return n_grid
################# Removing the SPORADIC grids ########################
def clean_grid(grid_list):
    print 'Cleaning the grid_list'
    beta = 0.3
    gkeys = grid_list.keys()
    for key in list(gkeys):
        value = grid_list[key]
        tg = value['tg']
        ct = t
        tm = value['tm']
        diff = ct-tg
        pi = (cl*(1-math.pow(delta,(diff+1))))/(N*(1-delta))
        if ((value['density'] < pi) and (t >= (1 + beta)*tm)):
            status = value['status']
            if status == 'SPORADIC' :
                if (t-tg >= gaptime) :
                    value['tg'] = 0
                    value['tm'] = t
                    value['density'] = 0
                    value['label']='NO_CLASS'   
                    value['status'] = "NORMAL"
            elif status == 'NORMAL' :
                value['status'] = "SPORADIC"
        else :
            value['status'] = "NORMAL"
        grid_list[key] = value       #i!

################ Adjusting the clusters ##############################
def adjust_cluster(grid_list,class_name):
    print "Adjusting the grids...."
    update_grids(grid_list,class_name)
    for g in grid_list:
        print 'Updation started'
        print "Grid:",g
        value = grid_list[g]
        tg = value['tg']
        prev_call = t- gaptime
        if prev_call < tg <=t:
            ckey = get_cluster(g)
            if ckey != None :
                c= cluster[ckey]
                if ( g_type(value) is "SPARSE" ):              ## Checking for SPARSE grid
                        print "Sparse grid"
                        if g not in c:
                            continue
                        else :
                            c.remove(g)                                             ## delete g from cluster c
                            value['label'] = 'NO_CLASS'
                            resolve_connectivity(ckey,class_name)
                elif ( g_type(value) is "DENSE" ) :                          ## checking for dense grid
                     print "Dense grid"
                     h = max_size_cluster(neigh_grid(g),cluster)
                     if h != None :
                         chkey = get_cluster(h)
                         if chkey != None :
                             hvalue = grid_list[h]
                             ch= cluster[chkey]
                             if (g_type(hvalue) is "DENSE") :
                                 if value['label'] == "NO_CLASS" :
                                     value['label'] = chkey
                                     cluster[chkey].append(g)
                                 else :
                                     g_cluster = get_cluster(g)
                                     if g_cluster == c:
                                         if len(c) > len(ch):
                                             move(chkey,ckey)
                                         elif len(c) <= len(ch):
                                             move(ckey,chkey)
                             elif (g_type(hvalue) is "TRANSITIONAL") :
                                 ctemp= list(ch)
                                 ctemp.append(g)
                                 
                                 if value['label'] == "NO_CLASS" and isoutside(h,ctemp):
                                     value['label'] = chkey    
                                     ch.append(g)
                                 elif len(c) >= len(ch):                   ##move: move h from cluster ch to c
                                     valueh= grid_list[h]
                                     valueh['label']= ckey
                                     c.append(h)
                                     if h not in ch:
                                         continue
                                     else :
                                         ch.remove(h)
                             cluster[chkey]=ch 
                             
                elif ( g_type(value) is "TRANSITIONAL" ) :
                    print "Transitional grid"
                    n_grid= neigh_grid(g)
                    h = max_size_cluster(n_grid,cluster)
                    if h != None :
                        ch= get_cluster(h)
                        if ch != None :
                            ctemp= list(cluster[ch])
                            while(n_grid):
                                n_grid.remove(h)
                                if isoutside(g,ctemp):
                                    #value= grid_list[g]
                                    l = value['label']
                                    value['label']= ch
                                    if l != 'NO_CLASS':
                                        c.remove(g)
                                        cluster[ch].append(g)                        
                                    break
                                h = max_size_cluster(n_grid,cluster)
                                if h != None :
                                    ch= get_cluster(h)
                                    if ch != None :
                                        ctemp= list(cluster[ch])
                                        ctemp.append(g)                                                        ## find largest c' satisfying that g is an outside grid among neighbouring clusters of g
                                    
                cluster[ckey]=c #i!
        grid_list[g]= value #i!

    #plot_cluster()
    return class_name
########## Checking for unconnected clusters and resolving it ##################         
def resolve_connectivity(ckey,class_name):
    print "Checking for unconnected clusters...."
    c2=list()
    for grid in cluster[ckey] :                                                         #checking unconnection
        value = grid_list[grid]
        if g_type(value) is "SPARSE":
            c2 = cluster[ckey]
            c2.remove(grid)
            cluster[ckey]=c2 #i!
            value['label'] = 'NO_CLASS'
            grid_list[grid]=value #i!
            break
        if isoutside(grid,cluster[ckey]) is False :
            break    
        else:
            for ngrid in neigh_grid(grid):
                if ngrid in cluster[ckey]:
                    break
                else:
                    value=grid_list[ngrid]
                    if g_type(value) is "DENSE" :
                        cluster[class_name]=[]
                        cluster[class_name].append(ngrid)
                        print 'Cluster', class_name,'created'
                        class_name += 1
############## Updating the density of grids in before adjust clustering #####################
def update_grids(grid_list,class_name) :
    for grid in grid_list :
        vector = grid_list[grid]
        update_charvec(grid,0)
        if(g_type(vector) == "DENSE"):
            l = vector.setdefault('label',class_name)
            if l == class_name :
                cluster[class_name]= []
                cluster[class_name].append(grid)
                print 'Cluster',class_name, 'created'
                class_name +=1
        else:
            l= vector.setdefault('label','NO_CLASS')
        grid_list[grid]=vector #i!
#################### Getting the cluster in which the grid exists ###################
def get_cluster(g) :
    key = grid_list[g].setdefault('label','NO_CLASS')
    if key!= 'NO_CLASS':
        return key
    else :
        return None
################## Getting the maximum size grid ##################
def max_size_cluster(n_grid,cluster):
    max_len = 0
    max_cluster_grid = None
    for key in n_grid :
        ch = get_cluster(key)
        if ch == None :
            print "Key has no cluster yet"
        else:
            size = len(cluster[ch])
            if size > max_len :
                max_len = size
    return max_cluster_grid
#############################################################################
############### Initializing the data structures and variables ##############
grid_list ={}
tgrid_list = {}
cluster = {}
input_list = []
class_name=1
start_time = datetime.datetime.now()
########### Reading the records from file to a list #################
input_list=pd.read_csv('test1.csv')  ## TODO: Normalize the data streams    
#######


cm=38
cl=0.67
p=50
delta=0.998

delta,dm,dl,gaptime,N,p=param_gen(delta,cm,cl,p)  #@#
grid_size =input_list.apply(lambda x:x.max()/p) #@# ##TODO: put 1 instead of x.max() if the data is normalized
tc= 0
t =1
for i in range(input_list.shape[0]):
    rec=input_list.iloc[i,:]
    g= get_den_grid(rec)  ## g stores the key for the input record    
    d =1000              ## TODO:constant factor to be multiplied to density 
    update_charvec(g,d) ##updates the density for the new input point
    timediff= t-tc
    if timediff == gaptime :
        print 'calling initial clustering'
        class_name = initial_clustering(grid_list,class_name)
        print 'cluster after initialisaiton', cluster
    if (timediff % gaptime)== 0:
        clean_grid(grid_list)
        class_name = adjust_cluster(grid_list,class_name)
        print 'after adjust cluster', cluster
    t +=1
