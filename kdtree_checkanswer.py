#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 15:44:06 2016

@author: root
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

def distsep(point,data):
    #the [ra,dec] of all stars should be in unit of degree
    #return angular seperation in unit of degree
    xp=point[1];yp=point[2];xd=data[1];yd=data[2]
    dist=np.sqrt((xd-xp)**2+(yd-yp)**2)

    return dist

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

def findneighbor(point,data,maxdtheta,mode):
    #point is [ra,dec] of a single star; data is position dataset
    #mode==1:distance seperation
    #mode==2:angular seperation
    if mode==1:
        dtheta=distsep(point,data)
    if mode==2:
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

def query(point,tree,maxtheta,mode):
    #suppose point is from cm1, kdtree is KDTREE generated through cm2
    #return [cm2_indx,dtheta]
    #in the official code maxtheta=2 arcsec=2/3600 degree
    if type(tree)==type(point):
        indx,dtheta=findneighbor(point,tree,maxtheta,mode)
        return indx,dtheta
    else:
        if point[tree.divdim]<tree.split:
            #print 'l'
            return query(point,tree.l,maxtheta,mode)           
        else:
            #print 'r'
            return query(point,tree.r,maxtheta,mode)
    
###############################################################################
cm1=np.array([m1,m2])#gaia result
cm2=np.array([match1,match2])#self result
maxtheta=0.0 #2 arcsec

data1=prepare(cm1);data2=prepare(cm2)
tree=kdtree(data1)#use gaia result to build tree
nstar=data2.shape[1]#use self result to inquery
mismatch_self=np.arange(0)
for counter in np.arange(nstar):
    point=pick(data2,counter)
    indx,dtheta=query(point,tree,maxtheta,1)
    if indx=='no neighbor':
        mismatch_self=np.append(mismatch_self,counter)
    #print indx,counter,dtheta

tree=kdtree(data2)#use self result to build tree
nstar=data1.shape[1]#use gaia result to inquery
mismatch_gaia=np.arange(0)
for counter in np.arange(nstar):
    point=pick(data1,counter)
    indx,dtheta=query(point,tree,maxtheta,1)
    if indx=='no neighbor':
        mismatch_gaia=np.append(mismatch_gaia,counter)
    #print indx,counter,dtheta

#Further checking the mismatch data in apogee catalog
mindx_gaia=m1[mismatch_gaia]
mindx_self=match1[mismatch_self]
ra_gaia=apogee_cat['RA'][mindx_gaia];dec_gaia=apogee_cat['DEC'][mindx_gaia]
ra_self=apogee_cat['RA'][mindx_self];dec_self=apogee_cat['DEC'][mindx_self]

cm1=np.array([ra_gaia,dec_gaia])#gaia result
cm2=np.array([ra_self,dec_self])#self result
maxtheta=0.0 #2 arcsec

data1=prepare(cm1);data2=prepare(cm2)
tree=kdtree(data1)#use gaia result to build tree
nstar=data2.shape[1]#use self result to inquery
mismatch_apogee=np.arange(0)
for counter in np.arange(nstar):
    point=pick(data2,counter)
    indx,dtheta=query(point,tree,maxtheta,2)
    if indx=='no neighbor':
        mismatch_apogee=np.append(mismatch_apogee,counter)
    #print indx,counter,dtheta   

#Further checking the mismatch data in tgas catalog
mindx_gaia=m2[mismatch_gaia]
mindx_self=match2[mismatch_self]
ra_gaia=tgas_cat['ra'][mindx_gaia];dec_gaia=tgas_cat['dec'][mindx_gaia]
ra_self=tgas_cat['ra'][mindx_self];dec_self=tgas_cat['dec'][mindx_self]

cm1=np.array([ra_gaia,dec_gaia])#gaia result
cm2=np.array([ra_self,dec_self])#self result
maxtheta=0.0 #2 arcsec

data1=prepare(cm1);data2=prepare(cm2)
tree=kdtree(data1)#use gaia result to build tree
nstar=data2.shape[1]#use self result to inquery
mismatch_tgas=np.arange(0)
for counter in np.arange(nstar):
    point=pick(data2,counter)
    indx,dtheta=query(point,tree,maxtheta,2)
    if indx=='no neighbor':
        mismatch_tgas=np.append(mismatch_tgas,counter)
    #print indx,counter,dtheta   

    
