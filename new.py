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
input_l = []
class_name = 1

################## Functions ###################################
############### Taking the parameters as input_l from user #################
def param_gen(delta,cm,cl,p):
    
    #global input_l
    
    ### EH: If N is too large show error and ask to decrease value of p. 
    #input_l = np.array(input_l)
    N=p**input_l.shape[1] ##shape[1] gives the dimension of the input_l i.e. the number of columns
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
    dm = 100*cm/(N*(1 - delta))
    dl = 1000*cl/(N*(1 - delta))
    gaptime = math.floor(d)  # TODO: Make sure gaptime  doesn't become zero # gaptime becomes 0 when cm is small and cl is large. #  d3 should be less than lambda
    #gaptime=20
    assert gaptime != 0, "Modulo by zero (timediff % gaptime) !!!"
    print 'Gaptime :', gaptime
    return delta,dm,dl,gaptime,N,p
    
########### Reading the records from file to a list #################     PROBLEM IN FUNCTION
def get_den_grid(rec):
    grid_size= [400/50,240000/50]
    temp=rec
    key=(temp/grid_size).apply(int)
    #key=(int((int)temp[0]/grid_size[0]),int((int)temp[1]/grid_size[1]))
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
        print "Moving grid",key,"from cluster",c1,"to cluster",c 
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
                print cluster
                class_name +=1
        else :
            grid_list[g]['label']= 'NO_CLASS'

    for c in list(cluster):          ## EH: list so as to avoid runtime error: dict changed size during iteration.
        if(cluster.get(c,0)):           ## 
            for g in cluster[c]:
                if(isOutside(g,c)):
                    for h in neigh_grid(g):
                        hvector= grid_list[h]
                        ch= hvector['label']
                        if (ch!= "NO_CLASS"):
                            if(ch!=c):      #EH: To avoid infinite loop where grid keeps moving in the same cluster.
                                if(len(cluster[c])> len(cluster[ch])):
                                    move(ch,c)
                                else :
                                    move(c,ch)
                        elif (g_type(hvector) == "TRANSITIONAL" ) :
                            print "Moving grid",h,"to cluster", c
                            grid_list[h]['label']= c
                            cluster[c].append(h)
                            print "Now class of", h,"is",grid_list[h]['label']
    return
    
################# Removing the SPORADIC grids ########################
def clean_grid():
    global grid_list
    print 'Cleaning the grid_list...'
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
                    print "Cleaning grid ",key
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
    print "Finished Cleaning..."
    
############## Updating the density of grids before adjusting clusters #####################
def update_grids() :
    global grid_list
    global class_name
    for grid in grid_list :
        update_charvec(grid,0)
        if(g_type(grid_list[grid]) == "DENSE"):
            l = grid_list[grid].setdefault('label',class_name)
            if l == class_name :
                cluster[class_name]= []
                cluster[class_name].append(grid)
                print 'Cluster',class_name, 'created'
                class_name+=1
        else:
            grid_list[grid]['label'] = 'NO_CLASS'
            
################## Getting the grid whose cluster has the maximum size ##################
def max_size_cluster_grid(n_grid):
    global cluster
    max_len = 0
    max_cluster_grid = None
    for key in n_grid :
        ch = get_cluster(key)
        if ch == None :
            print " "
        else:
            size = len(cluster[ch])
            if size > max_len :
                max_len = size
                max_cluster_grid=key
    return max_cluster_grid
    
#################### Getting the cluster in which the grid exists ###################
def get_cluster(g) :
    key = grid_list[g].setdefault('label','NO_CLASS')
    if key!= 'NO_CLASS':
        return key
    else :
        return None
        
########## Checking for unconnected clusters and resolving them into separate clusters ##################         
def resolve_connectivity(ckey):
    global class_name
    print "Checking for unconnected clusters...."
    
    #########################-------------TO BE WRITTEN -----------------################
#    c2=list()
#    for grid in cluster[ckey] :                                                         #checking unconnection
#        value = grid_list[grid]
#        if g_type(value) is "SPARSE":
#            c2 = cluster[ckey]
#            c2.remove(grid)
#            cluster[ckey]=c2 #i!
#            value['label'] = 'NO_CLASS'
#            grid_list[grid]=value #i!
#            break
#        if isoutside(grid,cluster[ckey]) is False :
#            break    
#        else:
#            for ngrid in neigh_grid(grid):
#                if ngrid in cluster[ckey]:
#                    break
#                else:
#                    value=grid_list[ngrid]
#                    if g_type(value) is "DENSE" :
#                        cluster[class_name]=[]
#                        cluster[class_name].append(ngrid)
#                        print 'Cluster', class_name,'created'
#                        class_name += 1
    
################ Adjusting the clusters ##############################
def adjust_cluster():
    global class_name
    global grid_list
    global cluster
    print "Adjusting the grids into appropriate clusters...."
    update_grids()
    
    for g in grid_list:
        print "Updation for grid ", g
        vector = grid_list[g]
        if(g_type(vector) is "SPARSE"):
            print "Sparse grid"
            ckey= get_cluster(g)
            if ckey!= None :
                c = cluster[ckey]
                if g in c:
                    cluster[ckey].remove(g)
                    grid_list[g]['label'] = 'NO_CLASS'
                    resolve_connectivity(ckey)
        elif(g_type(vector) is "DENSE"):
            print "Dense Grid"
            maxg = max_size_cluster_grid(neigh_grid(g))
            if maxg != None :
                chkey = get_cluster(maxg)
                if chkey != None :
                    if(g_type(grid_list[maxg]) is "DENSE"):
                        if vector['label'] == 'NO_CLASS' :
                            grid_list[g]['label'] = chkey
                            cluster[chkey].append(g)
                        else:
                            ckey = get_cluster(g)
                            if ckey != chkey :      ##EH : Only merge clusters when we have 2 different clusters.
                                if len(cluster[ckey]) > len(cluster[chkey]):
                                    move(chkey, ckey)
                                elif len(cluster[ckey]) <= len(cluster[chkey]) :
                                    move(ckey,chkey)                            
                    elif(g_type(grid_list[maxg]) is "TRANSITIONAL"):
                        if vector['label'] == 'NO_CLASS':
                            cluster[chkey].append(g)
                            grid_list[g]['label'] = chkey 
                            if isOutside(g,chkey) :
                                print " "
                            else :
                                grid_list[g]['label'] = 'NO_CLASS'
                                cluster[chkey].remove(g)
                        elif len(cluster[ckey]) >= len(cluster[chkey]):
                            if ckey != chkey :
                                grid_list[maxg]['label'] = ckey
                                cluster[ckey].append(maxg)
                                if maxg in cluster[chkey] :
                                    cluster[chkey].remove(maxg)
        elif(g_type(vector) is "TRANSITIONAL"):
            print "Transitional Grid"
            n_grid = neigh_grid(g)
            while(n_grid):
                maxg= max_size_cluster_grid(n_grid)
                if maxg != None :
                    chkey= get_cluster(maxg)
                    if chkey != None:
                        if isOutside(g,chkey):
                            if(vector['label'] != 'NO_CLASS'):
                                cluster[get_cluster(g)].remove(g)
                            cluster[chkey].append(g)
                            grid_list[g]['label'] = chkey
                            break
                    n_grid.remove(maxg)
    return                    
                        

    
########### Reading the records from file to a list #################
#text_file = open("test3.csv", "r")
#reader = csv.reader(text_file,delimiter = ',')
#for row in reader :
#    input_l.append(row)  ## TODO: Normalize the data streams
#start_time = datetime.datetime.now()   


input_l=pd.read_csv('test2.csv')


####################input_l parameters from the user###################
cm=55
cl=0.99
p=50
delta=0.998
maxima = lambda x:x.max()/p
grid_size = input_l.apply(maxima)
delta,dm,dl,gaptime,N,p=param_gen(delta,cm,cl,p)
 #@# ##TODO: put 1 instead of x.max() if the data is normalized

tc=0
t=1
for i in range(0,25): #(0,input_l_l.shape[0]):
    rec=input_l.iloc[i,:]   ## iloc gets the input_l row at location i 
    g= get_den_grid(rec)  ## g stores the key for the input_l_l record    
    d =1000              ## TODO:constant factor to be multiplied to density 
    update_charvec(g,d) ##updates the density for the new input_l_l point
    timediff= t-tc
    print "Processed record" , rec , "into grid ", g
    if timediff == gaptime :
        print 'Cluster before initial clustering', cluster
        initial_clustering()
        print 'Cluster after initial clustering', cluster
    if (timediff % gaptime)== 0:
        clean_grid()
        adjust_cluster()
        print 'Cluster after adjusting', cluster
    t +=1
