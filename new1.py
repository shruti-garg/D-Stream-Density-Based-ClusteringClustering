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
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib  import cm
#############################################################################
############### Initializing the data structures and variables ##############
grid_list ={}                   ##{} is for dictionary, [] is for list
cluster = {}
input_l = []
class_name = 1
new_class=1

################## Functions ###################################
############### Taking the parameters as input_l from user #################
def param_gen(delta,cm,cl,p):
    ### EH: If N is too large show error and ask to decrease value of p. 
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
    dm = 1000*cm/(N*(1 - delta))
    dl = 1000*cl/(N*(1 - delta))
    gaptime = math.floor(d)  # TODO: Make sure gaptime  doesn't become zero # gaptime becomes 0 when cm is small and cl is large. #  d3 should be less than lambda
    assert gaptime != 0, "Modulo by zero (timediff % gaptime) !!!"
    print 'Gaptime :', gaptime
    return delta,dm,dl,gaptime,N,p
########### Reading the records from file to a list #################     PROBLEM IN FUNCTION
def get_den_grid(rec):
    global grid_list
    temp=rec
    key=(temp/grid_size).apply(int)
    key=tuple(key)
    if key not in grid_list :
        grid_list[key]= {'tg':0,'tm':0, 'density':0,'status':'NORMAL','label':'NO_CLASS'}                                              ##   adds key to the grid with empty values
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
    elif (( D < dm ) and ( D > dl )) :
        gtype = "TRANSITIONAL"
    return gtype
############## Checking if the grid is an outside grid ##################### 
def isOutside(g,c):
    global cluster
    for ng in neigh_grid_all(g):        #consider all possible neigh grids to check if it's an outside grid.
        if ng not in cluster[c]:
            return True
    return False
############## Finding the neighbour grids which are in the grid_list ######################    
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
############## Finding the neighbour grids for all grids possible######################    
def neigh_grid_all(g) :
    n_grid = []
    for j in range(0,len(g)):
        temp1= list(g)
        temp2 = list(g)
        temp1[j] = temp1[j]-1
        temp1 = tuple(temp1)
        temp2[j] = temp2[j]+1
        temp2 = tuple(temp2)
        n_grid.append(temp1)
        n_grid.append(temp2)
    return n_grid
############# Move grids from one cluster to another #####################
def move(c1, c):
    global cluster
    global grid_list
    for key in cluster[c1]:
        print "Moving grid",key,"from cluster",c1,"to cluster",c 
        grid_list[key]['label'] = c
        if key not in cluster[c]:
            cluster[c].append(key)
    del cluster[c1]
    print "Cluster ", c1 , "deleted"
    return
################ Forming the initial clusters ###########################
def initial_clustering():
    global new_class
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
                for i in range(1,class_name+1):
                    if(i not in cluster):
                        new_class =i
                        break
                cluster[new_class]= []
                cluster[new_class].append(g)
                print 'Cluster',new_class, 'created'
                print cluster
                if new_class == class_name:
                    class_name +=1
        else :
            grid_list[g]['label']= 'NO_CLASS'

    for c in list(cluster):          ## EH: list so as to avoid runtime error: dict changed size during iteration.
        if(cluster.get(c,0)):           ## get returns value for given key, if key not present returns 0.
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
                                    break##EH: C cluster finishes, so no pint in iterating in inner for loop
                        elif (g_type(hvector) == "TRANSITIONAL" ) :
                            print "Moving grid",h,"to cluster", c
                            grid_list[h]['label']= c
                            if h not in cluster[c]:     ##EH: so that a cluster doesn't contain multiple grids with same key
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
                    if(value['label'] != 'NO_CLASS'):
                        cluster[get_cluster(key)].remove(key)
                        if(len(cluster[get_cluster(key)]) == 0):
                            del cluster[get_cluster(key)]
                            print "Deleting Cluster", get_cluster(key)
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
#        if(g_type(grid_list[grid]) == "DENSE"):
#            l = grid_list[grid].setdefault('label',class_name)
#            if l == class_name :
#                cluster[class_name]= []
#                cluster[class_name].append(grid)
#                print 'Cluster',class_name, 'created'
#                class_name+=1
#        else:
#            if(grid_list[grid]['label'] != 'NO_CLASS'):
#                cluster[get_cluster(grid)].remove(grid)
#                if(len(cluster[get_cluster(grid)]) == 0):
#                    del cluster[get_cluster(grid)]
#                    print "Deleting Cluster", get_cluster(grid)
#                grid_list[grid]['label'] = 'NO_CLASS'
                
################## Getting the grid whose cluster has the maximum size ##################
def max_size_cluster_grid(n_grid):
    global cluster
    max_len = 0
    max_cluster_grid = None
    for key in n_grid :
        ch = get_cluster(key)
        if ch == None :
            print " "
        elif ch in cluster:
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
def resolve_connectivity(ckey,g):
    global class_name
    global cluster
    global grid_list
    global new_class
    print "Checking for unconnected clusters...."
    #If g is an inside grid for cluster, do nothing.
    if isOutside(g, ckey) is False :
        return
    #for all grids in cluster, make new cluster for dense grids, remove sparse grids, label as NO_CLASS for transitional grids.
    for grid in list(cluster[ckey]) : 
        vector= grid_list[grid]
        if g_type(vector) is "SPARSE":
            if grid in cluster[ckey]:
                cluster[ckey].remove(grid)
                if(len(cluster[ckey]) == 0):
                    del cluster[ckey]
                print "Deleting Cluster", ckey
            grid_list[grid]['label']= 'NO_CLASS'
        elif g_type(vector) is "DENSE":
            for i in range(1,class_name+1):
                if(i not in cluster):
                    new_class =i
                    break
            cluster[new_class]= []
            cluster[new_class].append(grid)
            if(grid_list[grid]['label'] != 'NO_CLASS'):
                cluster[get_cluster(grid)].remove(grid)
                if(len(cluster[get_cluster(grid)]) == 0):
                    del cluster[get_cluster(grid)]
                    print "Deleting Cluster", get_cluster(grid)
            grid_list[grid]['label']= new_class
            print 'Cluster',new_class, 'created'
            if new_class == class_name:
                class_name +=1
        
        elif g_type(vector) is "TRANSITIONAL":
            grid_list[grid]['label'] = 'NO_CLASS'
    if ckey in cluster:
        del cluster[ckey]                                   # Delete the whole cluster, it would be adjusted accordingly.
    print "Deleting Cluster ",ckey
    return
################ Adjusting the clusters ##############################
def adjust_cluster():
    global class_name
    global new_class
    global grid_list
    global cluster
    update_grids()
    for g in grid_list:
        #print "Updation for grid ",g, 
        ckey = get_cluster(g)
        vector = grid_list[g]
        if(g_type(vector) is "SPARSE"):
            print "Sparse Grid", g
            ckey= get_cluster(g)
            if ckey!= None :
                c = cluster[ckey]
                if g in c:
                    cluster[ckey].remove(g)
                    if(len(cluster[ckey]) == 0):
                        del cluster[ckey]
                        print "Deleting Cluster", ckey
                    grid_list[g]['label'] = 'NO_CLASS'
                    resolve_connectivity(ckey,g)
        elif(g_type(grid_list[g]) is "DENSE"):
            print "Dense Grid", g
            maxg = max_size_cluster_grid(neigh_grid(g))
            if maxg != None :
                chkey = get_cluster(maxg)
                if chkey != None :
                    if(g_type(grid_list[maxg]) is "DENSE"):
                        if grid_list[g]['label'] == 'NO_CLASS' :
                            grid_list[g]['label'] = chkey
                            if g not in cluster[chkey]:             ##EH
                                cluster[chkey].append(g)
                        else:
                            ckey = get_cluster(g)
                            if ckey != chkey :      ##EH : Only merge clusters when we have 2 different clusters.
                                if len(cluster[ckey]) > len(cluster[chkey]):
                                    move(chkey, ckey)
                                elif len(cluster[ckey]) <= len(cluster[chkey]) :
                                    move(ckey,chkey)                            
                    elif(g_type(grid_list[maxg]) is "TRANSITIONAL"):
                        if grid_list[g]['label'] == 'NO_CLASS':
                            if g not in cluster[chkey]:
                                cluster[chkey].append(g)
                            grid_list[g]['label'] = chkey 
                            if isOutside(g,chkey) :
                                print " "
                            else :
                                grid_list[g]['label'] = 'NO_CLASS'
                                cluster[chkey].remove(g)
                                if(len(cluster[chkey]) == 0):
                                    del cluster[chkey]
                                    print "Deleting Cluster", chkey
                        elif (ckey in cluster and chkey in cluster) and (len(cluster[ckey]) >= len(cluster[chkey])): ##EH: execute only when ckey and chkey exist in cluster
                            if ckey != chkey :
                                grid_list[maxg]['label'] = ckey
                                if maxg not in cluster[ckey]:
                                    cluster[ckey].append(maxg)
                                if maxg in cluster[chkey] :
                                    cluster[chkey].remove(maxg)
                                    if(len(cluster[chkey]) == 0):
                                        del cluster[chkey]
                                        print "Deleting Cluster", chkey
            if(grid_list[g]['label'] == 'NO_CLASS'):                 #Correction: Make a new cluster with this dense grid if it cannot be merged into other cluster.
                for i in range(1,class_name+1):
                    if(i not in cluster):
                        new_class =i
                        break
                cluster[new_class]= []
                cluster[new_class].append(g)
                print 'Cluster',new_class, 'created'
                grid_list[g]['label']= new_class
                if new_class == class_name:
                    class_name +=1
        elif(g_type(grid_list[g]) is "TRANSITIONAL"):
            print "Transitional grid", g
            if(grid_list[g]['label'] != 'NO_CLASS'):
                if get_cluster(g) in cluster:
                    if g in cluster[get_cluster(g)]:
                        cluster[get_cluster(g)].remove(g)
                    if(len(cluster[get_cluster(g)]) == 0):
                        del cluster[get_cluster(g)]
                        print "Deleting Cluster", get_cluster(g)
                    grid_list[g]['label'] = 'NO_CLASS'
            n_grid = neigh_grid(g)
            while(len(n_grid)):
                maxg= max_size_cluster_grid(n_grid)
                if maxg == None:                            ## EH: To avoid infinite loop when no neighbourin ggrid belong to a cluster
                    break
                if maxg != None :
                    chkey= get_cluster(maxg)
                    if chkey != None:
                        if isOutside(g,chkey):
                            if g not in cluster[chkey]:
                                cluster[chkey].append(g)
                            grid_list[g]['label'] = chkey
                            #print "I'm changing label to ",chkey
                            break
                    n_grid.remove(maxg)
    return
##############################map every point to its cluster label#####################
#def map_label(x):
#    key=tuple((x/grid_size).apply(int))
#    clus=grid_list[key]['label']
#    return clus
#def get_label():
#    clus_l=[]
#    for i in range(0,input_l.shape[0]):
#        c=map_label(input_l.iloc[i,:])
#        clus_l.append(c)
#    return clus_l
def plot_it():
    #labels_series=pd.Series(get_label())
    if (input_l.shape[1]==2):
        dim2_plot()
    if (input_l.shape[1]==1):
        dim1_plot()
    if (input_l.shape[1]>2):
        pca_plot()
    return
def dim2_plot(labels_series):
    d=pd.DataFrame(grid_list)
    d=d.T.reset_index()
    d['density']=d.density.apply(float)
    d['label']=d.label.apply(lambda x:-1 if x=='NO_CLASS' else x).apply(int)
    # Cluster Plot
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.set_title("Cluster plot",fontsize=14)
    ax.set_xlabel("X",fontsize=12)
    ax.set_ylabel("Y",fontsize=12)
    ax.scatter(x=d.level_0,y=d.level_1,c=d.label, marker = 'o')
    #density plot
    d.plot.scatter(x='level_0',y='level_1',c='density',title='Density plot')
    return
def pca_plot():
    d=pd.DataFrame(grid_list)
    key_df=pd.DataFrame(list(d.T.index))
    pca=PCA(n_components=2)
    key_df_pca=pd.DataFrame(data=(pca.fit_transform(key_df)),columns=['PC1','PC2'])
    key_df_pca=pd.concat([(d.T)[['density','label']].reset_index(drop=True),key_df_pca],axis=1)
    key_df_pca['density']=key_df_pca.density.apply(float)
    key_df_pca['label']=key_df_pca.label.apply(lambda x:-1 if x=='NO_CLASS' else x).apply(int)
    # Cluster Plot
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.set_title("PCA Cluster plot",fontsize=14)
    ax.set_xlabel("PC1",fontsize=12)
    ax.set_ylabel("PC2",fontsize=12)
    ax.scatter(x=key_df_pca.PC1.values,y=key_df_pca.PC2.values,c=key_df_pca.label.values, marker = 'o', cmap = cm.jet )
    #density plot
    key_df_pca.plot.scatter(x='PC1',y='PC2',c='density',title='Density PCA plot')
    return

input_l=pd.read_csv('dim4.csv')
####################input_l parameters from the user###################
cm=30
cl=0.1
p=10
delta=0.997
maxima = lambda x:x.max()/p
grid_size = input_l.apply(maxima)                   ##TODO: put 1 instead of x.max() if the data is normalized
delta,dm,dl,gaptime,N,p=param_gen(delta,cm,cl,p)
tc=0
t=1
for i in range(input_l.shape[0]): #range(input_l.shape[0]):
    rec=input_l.iloc[i,:]   ## iloc gets the input_l row at location i 
    g= get_den_grid(rec)  ## g stores the key for the input_l_l record    
    d =1000              ## TODO:constant factor to be multiplied to density 
    update_charvec(g,d) ##updates the density for the new input_l_l point
    timediff= t-tc
    if timediff == gaptime :
        print 'Cluster before initial clustering', cluster
        initial_clustering()
        print 'Cluster after initial clustering', cluster
    if (timediff % gaptime)== 0:
        clean_grid()
        adjust_cluster()
        print '****************Cluster at t=',t, cluster
    t +=1
