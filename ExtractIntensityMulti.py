# -*- coding: utf-8 -*-
"""
Created on Sun Oct 22 16:37:22 2017

@author: dzhaojie
"""
#import HelpFunctionsForCellTracking as HFCT
import os
import numpy as np
from skimage import io




def extract_intensity(num_of_field):
    F=io.imread('C:/Users/Desktop/AVG_flatfield-5x-590 nm LED.tif')
    F=F.astype(np.float64)
    mga=16844.556
    C=F/mga
    def cmask(index,radius,array):
          a,b = index
          nx,ny = array.shape
          y,x = np.ogrid[-a:nx-a,-b:ny-b]
          mask = x*x + y*y <= radius*radius

          return(sum(array[mask]))
      
        
    for field_of_view in num_of_field:
        #HFCT.clear_all
        try:
            images_name=os.listdir('C:/Users/data-'+str(field_of_view))
        except:
            continue

        xyl=np.load('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/xyl.npy')
        uniquexyl=np.unique(xyl[:,4])
        cell_cell_dis=np.zeros((xyl.shape[0],1))
        
        
        
        for ui in range (uniquexyl.shape[0]-1):
            num_of_cell=uniquexyl[ui].astype(np.int32)
            try:
                X=np.load('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/'+str(num_of_cell)+'cell/X.npy')
                Y=np.load('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/'+str(num_of_cell)+'cell/Y.npy')
            except:
                continue
            
            r=np.zeros((X.shape[0],1))
            g=np.zeros((X.shape[0],X.shape[1]))
            for i in range (X.shape[1]):
                for k in range (xyl.shape[0]):
                    cell_cell_dis[k]=np.sqrt(  np.square(X[0,i]-xyl[k,0]) +   np.square(Y[0,i]-xyl[k,1])  )
                diss=np.sort(cell_cell_dis,axis=0)
                if diss[1]>=8:
                    r[i]=4
                else:
                    r[i]=np.floor(diss[1]/2)
            for j in range (X.shape[0]):
                AO=io.imread(os.path.join('C:/Users/data-'+str(field_of_view),images_name[j]))
                A=AO.astype(np.float64)
                A=A/C
                for i in range (X.shape[1]):
                    g[j,i]=cmask((X[j,i],Y[j,i]),r[i],A)
            np.save('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/'+str(num_of_cell)+'cell/g.npy',g)