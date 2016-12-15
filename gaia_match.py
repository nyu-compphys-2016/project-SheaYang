#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 10:55:07 2016

@author: tao
"""

import gaia_tools.load as gload
from gaia_tools import xmatch
import matplotlib.pyplot as plt
import numpy as np

import time

tgas_cat=gload.tgas()
apogee_cat=gload.apogee()

tic = time.clock()
m1,m2,sep= xmatch.xmatch(apogee_cat,tgas_cat,colRA2='ra',colDec2='dec',epoch1=2000.,epoch2=2015.,swap=True)
toc = time.clock()
dt_astropy=toc - tic

plt.figure(1)
plt.xlabel('J')
plt.ylabel('delta_theta(arcsec)')
plt.scatter(apogee_cat['J'][m1],sep*3600,c='k',s=0.5,marker='.')

delta_ra_gaia=apogee_cat['RA'][m1]-tgas_cat['ra'][m2]
delta_dec_gaia=apogee_cat['DEC'][m1]-tgas_cat['dec'][m2]
plt.figure(2)
plt.xlabel('delta_dec (arcsec)')
plt.ylabel('delta_ra (arcsec)')
plt.scatter(delta_dec_gaia*3600,delta_ra_gaia*3600,c='k',s=0.5,marker='.')

