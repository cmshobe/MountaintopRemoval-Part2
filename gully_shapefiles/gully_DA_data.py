#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 30 14:13:38 2022

@author: sambower
"""

#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatch
import geopandas as gpd
import mplcursors
import seaborn as sbn
import matplotlib.image as img



path = '/Users/sambower/Documents/WVU/MTR/Data/clipped_by_watershed/gully_DA/'


bc = gpd.read_file(path+'bencreek_points.shp', delimiter = ',')
sf = gpd.read_file(path+'sprucefork_points.shp', delimiter = ',')
lc = gpd.read_file(path+'laurelcreek_points.shp', delimiter = ',')
wo = gpd.read_file(path+'whiteoak_points.shp', delimiter = ',')

bc.rename(columns = {'bencreekslo':'slope'}, inplace = True)
sf.rename(columns = {'spruceforks':'slope'}, inplace = True)
lc.rename(columns = {'laurelcreek':'slope'}, inplace = True)
wo.rename(columns = {'mudriverslo':'slope'}, inplace = True)

bc['WID'] = 'bencreek'
sf['WID'] = 'sprucefork'
lc['WID'] = 'laurelcreek'
wo['WID'] = 'whiteoak'

# get min and max precipitation values from NASA GCM data

path2 = '/Users/sambower/Documents/WVU/MTR/Data/climate/'

precip = pd.read_csv(path2+'yearly_data.csv', delimiter = ',')
precip = precip.transpose()
precip_array = precip.to_numpy()
precip_min = np.min(precip_array[2:,:]) #mm/yr
precip_max = np.max(precip_array[2:,:])

#%%

# DEFINE VARIABLES

#RANGE OF PRECIPITATION VALUES (P)

P_min = precip_min/1000 #m/yr
P_max = precip_max/1000 #m/yr

#RANGE OF EROSION VALUES EXTRACTED FROM GULLIES

E_min = 0.01 #m/yr
E_max = 1.54 #m/yr

#scaling factors

m = 0.5
n = 1

# min and max slope values 

slopes = pd.DataFrame()
slopes['slopes'] = slopes.append((sf['slope'],wo['slope'],bc['slope'],lc['slope']))
slopes['slopes'] = np.tan(np.deg2rad(slopes.slopes))

area = pd.DataFrame()
area['DA'] = area.append((sf['FLOW'],wo['FLOW'],bc['FLOW'],lc['FLOW']))

ids = pd.DataFrame()
ids['WID'] = ids.append((sf['WID'],wo['WID'],bc['WID'],lc['WID']))

df = pd.DataFrame()
df['WID'] = ids.WID
df['DA'] = area.DA
df['slope'] = slopes.slopes
df['product'] = df.DA * df.slope

max_point = df['product'].max()

S_min = df['slope'][df['product'] == df['product'].min()].to_numpy() #m/m
S_max = df['slope'][df['product'] == df['product'].max()].to_numpy() #m/m


A_min = df['DA'][df['product'] == df['product'].min()].to_numpy() #m^2
A_max = df['DA'][df['product'] == df['product'].max()].to_numpy() #m^2

#NOW do the min and max calculation for the highest and lowest possible erodibility scenarios

K_min = E_min/ ((P_max * A_max)**m * S_max**n)

K_max = E_max/ ((P_min * A_min)**m * S_min**n)


#%%

E_avg = 0.114 #m/yr (the average gully depth/average mine age (30yr)) #(E_min+E_max) /2
P_avg = np.mean(precip_array[2:,:])

#df['K_avg'] = abs(np.log10(E_avg/ ((P_avg * df['DA'])**m * df['slope']**n)))
df['K_avg'] = (E_avg/ ((P_avg * df['DA'])**m * df['slope']**n))

size = df['K_avg'].to_numpy()

df['color'] = None
df['color'][df['WID'] == 'bencreek'] = 'tab:blue'
df['color'][df['WID'] == 'sprucefork'] = 'gold'
df['color'][df['WID'] == 'laurelcreek'] = 'springgreen'
df['color'][df['WID'] == 'whiteoak'] = 'midnightblue'


fig, ax = plt.subplots(figsize = (10,10))
ax.grid(which='major', axis='both', zorder= 1)

ax.scatter(x = np.log10(df['DA'][df['WID'] == 'bencreek']),
           y = df['slope'][df['WID'] == 'bencreek'],
           s = (df['K_avg'][df['WID'] == 'bencreek'])*1000000,
           c = 'tab:blue',
           alpha = 0.50,
           edgecolor = 'k',
           marker = 's', zorder = 2)
ax.scatter(x = np.log10(df['DA'][df['WID'] == 'sprucefork']),
           y = df['slope'][df['WID'] == 'sprucefork'],
           s = (df['K_avg'][df['WID'] == 'sprucefork'])*1000000,
           c = 'gold',
           alpha = 0.50,
           edgecolor = 'k',
           marker = "^", zorder = 2)
ax.scatter(x = np.log10(df['DA'][df['WID'] == 'laurelcreek']),
           y = df['slope'][df['WID'] == 'laurelcreek'],
           s = (df['K_avg'][df['WID'] == 'laurelcreek'])*1000000,
           c = 'springgreen',
           alpha = 0.50,
           edgecolor = 'k',
           marker = '*', zorder = 2)
ax.scatter(x = np.log10(df['DA'][df['WID'] == 'whiteoak']),
           y = df['slope'][df['WID'] == 'whiteoak'],
           s = (df['K_avg'][df['WID'] == 'whiteoak'])*1000000,
           c = 'midnightblue',
           alpha = 0.50,
           edgecolor = 'k', zorder = 2)

ax.legend(['Ben Creek','Spruce Fork','Laurel Creek','White Oak'], fontsize = 16, loc = 'upper right')
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlabel('$\log_{10}drainage\hspace{.5}area\hspace{.5}(m^2)$',fontsize = 16)
ax.set_ylabel('slope $(m/m)$',fontsize = 16)

x_tail = 1.4
y_tail = 0.08
x_head = 1.4
y_head = 0.175
arrow = mpatch.FancyArrowPatch((x_tail, y_tail), (x_head, y_head),
                                  mutation_scale=20,
                                  alpha = .3,
                                  color = 'k',
                                  edgecolor = 'k')

ax.scatter(1.4,0.07,s= 20, c = 'k', zorder = 2)
ax.scatter(1.4,0.2,s= 700, c = 'k', zorder = 2)
ax.add_patch(arrow)
ax.text(1.2,0.00,'relative\nerodibility $K$', fontsize = 12, zorder = 2)
ax.text(1.5,0.065,'$2.5^{-5}$', zorder = 2)
ax.text(1.55,0.195,'$10^{-3}$', zorder = 2)
rec = mpatch.Rectangle((1.1,-0.01),0.9,0.25,fill = None, alpha = 0.25, edgecolor = 'k')
ax.add_patch(rec)

mplcursors.cursor()

#%%
a = df['K_avg'].to_numpy()
(n, bins, patches) = plt.hist(a, bins=20, label='hst')
mplcursors.cursor()
good_K_min = bins[0]
good_K_max = bins[3]

fig1, ax = plt.subplots(figsize = (10,10))
ax.grid(which='major', axis='both', zorder= -1.0)
ax = sbn.kdeplot(data = df['K_avg'], color = 'green', lw = 3, fill = True, zorder = 2)
# ax1 = ax.twinx()
# ax1.hist(df['K_avg'], bins = 20, zorder = 1)
ax.vlines(x = good_K_min, ymin = 0, ymax = 2500, ls = '--', lw = 2, color = 'k')
ax.vlines(x = good_K_max, ymin = 0, ymax = 2500, ls = '--', lw = 2, color = 'k')
ax.set_ylim(0,2500)
ax.set_xlabel('Possible K values', fontsize = 20)
ax.set_ylabel('Density', fontsize = 20)

ax.tick_params(axis='both', which='major', labelsize=16)
# ax.ticklabel_format(axis = 'x', style = 'sci', scilimits = (0,0), useMathText = True)
mplcursors.cursor()

#%%

K_min = good_K_min
K_max = good_K_max

start_time = 0
end_time = 200

coef = 0.0021

array = np.zeros(end_time)
array[0] = K_min

for i in range(1,end_time):
    array[i] = (coef)*array[i-1]**0.25
    # array[i] = 0.034*np.e**(-0.15*i)
    
k_array = K_max - array  

fig, ax = plt.subplots()
ax.scatter(0,good_K_max)
ax.scatter(30,good_K_min)
x = np.arange(0,end_time,1)
k_array = K_max - array
ax.plot(x,k_array)


time = 300
testshape = np.zeros(100)
testscape = np.zeros(len(testshape))
testshape[:50] = 0
testshape[50:] = 1

k_grids = np.zeros((len(testscape),time))
erosivity = np.zeros(time)
for n in range(time):
    if n < len(k_array):
        erosivity[n] = k_array[n]
        testscape[testshape == 0] = 9e-6
        testscape[testshape == 1] = k_array[n]
        k_grids[:,n] = testscape
        
    elif n < 200:
        testscape[testshape == 0] = 9e-6
        testscape[testshape == 1] = K_min
        k_grids[:,n] = testscape
        
    else:
        testscape[testshape == 0] = 9e-6
        testscape[testshape == 1] = 9e-6
        k_grids[:,n] = testscape


#%%


fig = plt.figure(constrained_layout=True, figsize=(12, 12))
spec = fig.add_gridspec(4, 4)

#PANEL (A)
ax1 = fig.add_subplot(spec[:2,:2])
ax1.set_xticks([])
ax1.set_yticks([])
im = img.imread('/Users/sambower/Desktop/gully.png')
ax1.axis('off')
ax1.imshow(im)
#THIS ONE IS FOR THE GULLY MAP --> LEAVE BLANK

#PANEL (B)
ax2 = fig.add_subplot(spec[:2,2:])
# ax2.grid(which='major', axis='both', zorder= 1)
ax2.scatter(x = np.log10(df['DA'][df['WID'] == 'bencreek']),
           y = df['slope'][df['WID'] == 'bencreek'],
           s = (df['K_avg'][df['WID'] == 'bencreek'])*2e5,
           c = 'tab:blue',
           alpha = 0.50,
           edgecolor = 'k',
           marker = 's', zorder = 2)
ax2.scatter(x = np.log10(df['DA'][df['WID'] == 'sprucefork']),
           y = df['slope'][df['WID'] == 'sprucefork'],
           s = (df['K_avg'][df['WID'] == 'sprucefork'])*2e5,
           c = 'gold',
           alpha = 0.50,
           edgecolor = 'k',
           marker = "^", zorder = 2)
ax2.scatter(x = np.log10(df['DA'][df['WID'] == 'laurelcreek']),
           y = df['slope'][df['WID'] == 'laurelcreek'],
           s = (df['K_avg'][df['WID'] == 'laurelcreek'])*2e5,
           c = 'k',
           alpha = 0.50,
           edgecolor = 'k',
           marker = '*', zorder = 2)
ax2.scatter(x = np.log10(df['DA'][df['WID'] == 'whiteoak']),
           y = df['slope'][df['WID'] == 'whiteoak'],
           s = (df['K_avg'][df['WID'] == 'whiteoak'])*2e5,
           c = 'midnightblue',
           alpha = 0.50,
           edgecolor = 'k', zorder = 2)
ax2.legend(['Ben Creek','Spruce Fork','Laurel Creek','White Oak'], fontsize = 10, loc = 'upper right')
ax2.tick_params(axis='both', which='major', labelsize=12)
ax2.set_xlabel('$\log_{10}drainage\hspace{.5}area\hspace{.5}(m^2)$',fontsize = 14)
ax2.set_ylabel('slope $(m/m)$',fontsize = 14)
x_tail,y_tail,x_head,y_head = 1.4,0.08,1.4,0.175
arrow = mpatch.FancyArrowPatch((x_tail, y_tail), (x_head, y_head),
                                  mutation_scale=20,
                                  alpha = .3,
                                  color = 'k',
                                  edgecolor = 'k')

ax2.scatter(1.4,0.07,s= 20, c = 'k', zorder = 2)
ax2.scatter(1.4,0.2,s= 200, c = 'k', zorder = 2)
ax2.add_patch(arrow)
ax2.text(1.2,0.00,'relative\nerodibility $K$', fontsize = 8, zorder = 2)
ax2.text(1.5,0.065,'$2.5^{-5}$', zorder = 2)
ax2.text(1.55,0.195,'$10^{-3}$', zorder = 2)
rec = mpatch.Rectangle((1.1,-0.01),0.9,0.25,fill = None, alpha = 0.25, edgecolor = 'k')
ax2.add_patch(rec)

#PANEL (C)
ax3 = fig.add_subplot(spec[2:, :2])
# ax3.grid(which='major', axis='both', zorder= -1.0)
ax3 = sbn.kdeplot(data = df['K_avg'], color = 'firebrick', lw = 3, fill = True, zorder = 2)
# ax1 = ax.twinx()
# ax1.hist(df['K_avg'], bins = 20, zorder = 1)
ax3.vlines(x = good_K_min, ymin = 0, ymax = 2500, ls = '--', lw = 2, color = 'k')
ax3.vlines(x = good_K_max, ymin = 0, ymax = 2500, ls = '--', lw = 2, color = 'k')
ax3.set_ylim(0,2500)
ax3.set_xlabel('Possible K values', fontsize = 14)
ax3.set_ylabel('Density', fontsize = 14)
ax3.tick_params(axis='both', which='major', labelsize=12)

#PANEL (D)
ax4 = fig.add_subplot(spec[2, 2:])
xd = np.arange(0,30, 0.25)
yd = xd**0.3
ax4.plot(xd,yd, lw = 3, color = 'mediumseagreen')
ax4.set_ylim(0,3.5)
ax4.set_xlim(0,30)
ax4.hlines(3,-5,35, ls = '--', color = 'k', lw = 2)
ax4.text(1,2.7, 'un-mined NDVI', color = 'k')
ax4.text(10,1.75, 'mined NDVI recovery ($t^{0.25}$)', color = 'mediumseagreen',rotation = 6)
ax4.set_yticks([])
ax4.set_ylabel('-   NDVI   +', fontsize = 14)
ax4.set_xlabel('time after reclamation (y)', fontsize = 14)
#PANEL (E)
ax5 = fig.add_subplot(spec[3, 2:])
# ax5.plot(erosivity)
xe = np.arange(0,300)
yd = xe**(-0.3)
yd[200:] = 0.1
ax5.plot(yd, lw = 3, zorder = 1)
ax5.set_yticks([])
ax5.set_xlim(0,300)
ax5.set_ylim(-.2,1.2)
ax5.set_ylabel('-  Erodibility $K$  +', fontsize = 14)
ax5.hlines(0.1,-5,350, ls = '--', lw = 2, color = 'k', zorder = 0)
ax5.hlines(1,-5,350, ls = '--', lw = 2, color = 'k', zorder = 0)
ax5.axvspan(0, 200, alpha=0.25, color='mediumseagreen', zorder = -1)
ax5.axvspan(200, 300, alpha=0.25, color='plum', zorder = -1)
ax5.text(20,-.05,'$K_{min}$')
ax5.text(20,1.05,'$K_{max}$')
ax5.text(100,0.45,'vegetation\nrecovery\nperiod', color = 'mediumseagreen')
ax5.text(225,0.45,'constant\nerodibility\nperiod', color = 'plum')
ax5.set_xlabel('time after reclamation (y)', fontsize = 14)
ax5.text(270,.15,'...10ky')
plt.savefig('/Users/sambower/Documents/WVU/MTR/Figures/K_figure.png', dpi = 300)









