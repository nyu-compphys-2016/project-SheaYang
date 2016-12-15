# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 22:03:31 2016

@author: sy182
"""
import gaia_tools.load as gload
import apogee.tools.read as apread
import matplotlib.pyplot as plt
import numpy as np

import time

import sys
sys.setrecursionlimit(10000)

# in this code, the original data should be [ra,dec]. prepare(data), all data be processed through other functions are [indx,ra,dec]
class node:
    def __init__(self):
        self.l=[]
        self.r=[]
        self.split='NAN'
        self.divdim='NAN'
    def set_left(self,data):
        self.l=data
    def set_right(self,data):
        self.r=data
    def set_split(self,split):
        self.split=split
    def set_divdim(self,divdim):
        self.divdim=divdim
        
def standard_split(data):
    #first column of data is indx, second column is ra, third is dec
    median=np.median(data,axis=1)
    dimlength=np.max(data,axis=1)-np.min(data,axis=1)
    indx=np.argmax(dimlength[1:])+1
    return median, indx
    
def midpoint_split(data):
    #this is only two dimensional data, so there will always be data on each side of the hyperplane
    #no need to slide the hyperplane
    #the tree is unbalanced but better constructed for inquiry
    median=(np.max(data,axis=1)+np.min(data,axis=1))/2
    dimlength=np.max(data,axis=1)-np.min(data,axis=1)
    indx=np.argmax(dimlength[1:])+1
    return median, indx
        

def prepare(data):
    #prepare [ra,dec] data pair (add index) before putting it in kdtree
    data=np.array(data)
    nstar=data.shape[1]
    data=np.array([np.arange(nstar),data[0],data[1]])   
    return data

def pick(data,indx):
    return np.array([data[0][indx],data[1][indx],data[2][indx]])

def angsep(point,data):
    #the [ra,dec] of all stars should be in unit of degree
    #return angular seperation in unit of degree
    rap=point[1]/180*np.pi;decp=point[2]/180*np.pi
    rad=data[1]/180*np.pi;decd=data[2]/180*np.pi
#    step1=np.cos(dec2)**2*np.sin(ra2-ra1)**2+(np.cos(dec1)*np.sin(dec2)-np.sin(dec1)*np.cos(dec2)*np.cos(ra2-ra1))**2
#    step2=np.sin(dec1)*np.sin(dec2)+np.cos(dec1)*np.cos(dec2)*np.cos(ra2-ra1)
#    step3=180/np.pi*np.arctan(np.sqrt(step1)/step2)
    step1=np.cos(decd)**2*np.sin(rad-rap)**2+(np.cos(decp)*np.sin(decd)-np.sin(decp)*np.cos(decd)*np.cos(rad-rap))**2
    step2=np.sin(decp)*np.sin(decd)+np.cos(decp)*np.cos(decd)*np.cos(rad-rap)
    step3=180/np.pi*np.arctan(np.sqrt(step1)/step2)
    return np.fabs(step3)

def findneighbor(point,data,maxdtheta):
    #point is [ra,dec] of a single star; data is position dataset
    dtheta=angsep(point,data)

    indx=np.argmin(dtheta)
    if dtheta[indx]>maxdtheta:
        return 'no neighbor','no neighbor'
    else:
        return data[0][indx], dtheta[indx]
    
def kdtree(data):
    if data.shape[1]<=10 :  
        return data
    else:
        branch=node()
        #median, indx_div=standard_split(data)
        median, indx_div=midpoint_split(data)
        mask_l=np.where(data[indx_div,:]<median[indx_div])
        mask_r=np.where(data[indx_div,:]>=median[indx_div])
        data_l=np.array([data[0][mask_l],data[1][mask_l],data[2][mask_l]])
        data_r=np.array([data[0][mask_r],data[1][mask_r],data[2][mask_r]])
        #print data_l.shape[1],data_r.shape[1]
        branch.set_left(kdtree(data_l))
        branch.set_right(kdtree(data_r))
        branch.set_split(median[indx_div])
        branch.set_divdim(indx_div)
        return branch

def query(point,tree,maxtheta):
    #suppose point is from cm1, kdtree is KDTREE generated through cm2
    #return [cm2_indx,dtheta]
    #in the official code maxtheta=2 arcsec=2/3600 degree
    if type(tree)==type(point):
        indx,dtheta=findneighbor(point,tree,maxtheta)
        return indx,dtheta
    else:
        if point[tree.divdim]<tree.split:
            #print 'l'
            return query(point,tree.l,maxtheta)           
        else:
            #print 'r'
            return query(point,tree.r,maxtheta)
    
###############################################################################
apogee_epoch=2000;tgas_epoch=2015
depoch=tgas_epoch-apogee_epoch
# Use proper motion to get both catalogs at the same time
dra=tgas_cat['pmra']/np.cos(tgas_cat['dec']/180.*np.pi)/3600000.*depoch
ddec= tgas_cat['pmdec']/3600000.*depoch

ra1=apogee_cat['RA'];dec1=apogee_cat['DEC']
ra2=tgas_cat['ra']-dra;dec2=tgas_cat['dec']-ddec

cm1=np.array([ra1,dec1])
cm2=np.array([ra2,dec2])
maxtheta=2.0/3600.0 #2 arcsec

tic = time.clock()
data1=prepare(cm1);data2=prepare(cm2)
tree=kdtree(data1)

nstar=data2.shape[1]
match1=np.arange(0);match2=np.arange(0);delta_theta=np.arange(0)
for counter in np.arange(nstar):
    point=pick(data2,counter)
    indx,dtheta=query(point,tree,maxtheta)
    if indx!='no neighbor':
        match1=np.append(match1,indx)
        match2=np.append(match2,counter)
        delta_theta=np.append(delta_theta,dtheta)
    #print indx,counter,dtheta
match1=match1.astype(int);match2=match2.astype(int)
toc = time.clock()
dt_self=toc - tic

delta_ra_self=apogee_cat['RA'][match1]-tgas_cat['ra'][match2]
delta_dec_self=apogee_cat['DEC'][match1]-tgas_cat['dec'][match2]
plt.figure(1)
plt.xlabel('delta_dec (arcsec)')
plt.ylabel('delta_ra (arcsec)')
plt.scatter(delta_dec_self*3600,delta_ra_self*3600,c='k',s=0.5,marker='.')

plt.figure(2)
plt.xlabel('J')
plt.ylabel('delta_theta (arcsec)')
plt.scatter(apogee_cat['J'][match1],delta_theta*3600,c='k',s=0.5,marker='.')


