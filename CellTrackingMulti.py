# -*- coding: utf-8 -*-
"""
@author: dzhaojie

"""
from __future__ import division
import os
import numpy as np
from skimage import io
from skimage import filters
import skimage
from skimage.morphology import extrema
import HelpFunctionsForCellTracking as HFCT


def cell_track(num_of_field):
    for field_of_view in num_of_field:
     
        
        HFCT.clear_all
        try:
            images_name=os.listdir('C:/Users/data-'+str(field_of_view))
        except:
            continue
        
        min_d=20 # minimum distance for single-cell and cell cluster tracking
        
        A=io.imread(os.path.join('C:/Users/data-'+str(field_of_view),images_name[0]))
        A=A.astype(np.float64)
        A=A-1200 # subtract certain amount of background noise
        A[A<0]=0
        NT=np.uint16(np.std(A))
        bw=extrema.h_maxima(A, NT)
        L=skimage.measure.label(bw)
        props=skimage.measure.regionprops(L,A)
        xyl=np.zeros((len(props),6))
        for i,properties in enumerate(props):
            xyl[i,1],xyl[i,0]=properties.centroid
            xyl[i,2]=properties.area
            xyl[i,3]=np.sum(properties.intensity_image)
     #use xyl[i,4] recode how many number of cells in the droplet that the cell i is in   
     #use xyl[i,5] recode the droplet number that the cell i is in
    
    
          
    ######################
    #Cell cluster recognization
    ######################
        A_cluster=io.imread(os.path.join('C:/Users/data-'+str(field_of_view),images_name[0]))
        B_cluster=io.imread(os.path.join('C:/Users/data-'+str(field_of_view),images_name[0+512])) #the number of images is denpent on each data set
        threshold=filters.threshold_otsu(B_cluster)
        bw_cluster=threshold<B_cluster
        L_cluster=skimage.measure.label(bw_cluster)
        props_cluster=skimage.measure.regionprops(L_cluster,A_cluster)
        xyl_cluster=np.zeros((len(props_cluster),6))
        for i,properties in enumerate(props_cluster):
            xyl_cluster[i,1],xyl_cluster[i,0]=properties.centroid
            xyl_cluster[i,2]=properties.area
            xyl_cluster[i,3]=np.sum(properties.intensity_image)
     #use xyl_cluster[i,5] recode the droplet number that the cell cluster i is in        
    
    
    
    ######################
    #Droplet recognization
    ######################        
        B2=io.imread(os.path.join('C:/Users/data-'+str(field_of_view),images_name[1024])) #the number of images is denpent on each data set
        threshold=filters.threshold_otsu(B2)
        bw2=threshold<B2
        L2=skimage.measure.label(bw2,neighbors=4) #use 4 connector otherwise serveral droplet will be recognized as one
        props2=skimage.measure.regionprops(L2,B2)
        xyl2=np.zeros((len(props2),5))
        for i,properties in enumerate(props2):
            xyl2[i,1],xyl2[i,0]=properties.centroid
            xyl2[i,2]=properties.equivalent_diameter
            xyl2[i,3]=properties.eccentricity
        row_indx=np.where(np.logical_or(xyl2[:,2]<30,xyl2[:,2]>500))
        xyl2=np.delete(xyl2,row_indx,axis=0)
        Droplet_Diameter=xyl2[:,2]
    #use xyl2[i,4] to recode how many cells this droplet has     
    
    
    
    ######################
    # Link the cells to the droplets
    ######################
        droplet_displace_X=0;#sometimes white field images is offset
        droplet_displace_Y=0;
        
        for i in range (xyl.shape[0]):
            if xyl[i,2] <= 20:
                for drop in range (xyl2.shape[0]):
                    if xyl2[drop,3]<=0.6:
                        drop_cell_dis=np.sqrt( np.square( xyl[i,0]-(xyl2[drop,0]+ droplet_displace_X) ) +np.square( xyl[i,1]-(xyl2[drop,1]+ droplet_displace_Y) ) )
                        if drop_cell_dis<=xyl2[drop,2]/2:
                            xyl[i,5]=drop
                            xyl2[drop, 4]=xyl2[drop,4]+1
                    else:
                        continue
            else:
                continue
                    
                    
                    
    ######################
    # Link the cell cluster to the droplets
    ######################    
        for i in range (xyl_cluster.shape[0]):     
            for drop in range (xyl2.shape[0]):
                if xyl2[drop,3]<=0.6:
                    drop_cluster_dis=np.sqrt(np.square( xyl_cluster[i,0]-(xyl2[drop,0]+ droplet_displace_X) ) +np.square( xyl_cluster[i,1]-(xyl2[drop,1]+ droplet_displace_Y) ) )
                    if drop_cluster_dis<=xyl2[drop,2]/2:
                        xyl_cluster[i,5]=drop
                else:
                    continue
        
       
        
    ######################
    # indentify single cells or multiple cells
    ######################
        for i in range (xyl.shape[0]):
            drop_label=xyl[i,5].astype(int)
            if drop_label==0:
                xyl[i,4]=900
            else:
                xyl[i,4]=xyl2[drop_label,4]
    
    
    
    ######################
    # save xyl and Droplet diameter
    ######################
        HFCT.make_path('C:/Users/data-results/local maxium separate cell/'+str(field_of_view))    
        np.save('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/xyl.npy',xyl)
        np.save('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/Droplet_Diameter.npy',Droplet_Diameter)
    
    
    
    ####################################################
    # start cell tracking, first is single-cell tracking
    ####################################################
        uniquexyl=np.unique(xyl[:,4])
        for ui in range (1):
            sort_cell=uniquexyl[ui]
            HFCT.make_path('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/'+str(ui+1)+'cell')
            N1=0
            j=0
            
            for Nt in range (xyl.shape[0]):
                if xyl[Nt,4]==sort_cell:
                    N1=N1+1
                else:
                    continue
            X=np.zeros((np.int32((len(images_name)-1)/2),N1))
            Y=np.zeros((np.int32((len(images_name)-1)/2),N1))
            a=np.zeros((np.int32((len(images_name)-1)/2),N1))
            g=np.zeros((np.int32((len(images_name)-1)/2),N1))
            DI=np.zeros((1,N1))
            N2=0
            for Nt in range (xyl.shape[0]):
                if xyl[Nt,4]==1:
                    N2=N2+1
                    X[j,N2-1] = xyl[Nt,0]
                    Y[j,N2-1] = xyl[Nt,1]
                    a[j,N2-1] = xyl[Nt,2]
                    g[j,N2-1] = xyl[Nt,3]
                    DI[j,N2-1] = xyl[Nt,5]
                else:
                    continue
                    
                    
            if N1>0:
                for j in range (1,np.int32((len(images_name)-1)/2)):
                        A=io.imread(os.path.join('C:/Users/data-'+str(field_of_view),images_name[j]))
                        A=A.astype(np.float64)
                        A=A-1200 # subtract certain amount of background noise
                        A[A<0]=0
                        NT=np.uint16(np.std(A))
                        if NT<=300:
                            NT=NT*2
                        bw=extrema.h_maxima(A, NT)
                        L=skimage.measure.label(bw)
                        props=skimage.measure.regionprops(L,A)
                        x=np.zeros((len(props),1))
                        y=np.zeros((len(props),1))
                        ma=np.zeros((len(props),1))
                        mg=np.zeros((len(props),1))
                        for i, properties in enumerate(props):
                            x[i]=properties.centroid[1]
                            y[i]=properties.centroid[0]
                            ma[i]=properties.area
                            mg[i]=np.sum(properties.intensity_image)
                        for m in range (N1):
                            dis=np.zeros((len(props),1))
                            for n in range (len(props)):
                                dis[n]=np.sqrt(  np.square(X[j-1,m]-x[n]) + np.square(Y[j-1,m]-y[n])   )
                            if np.min(dis)> min_d:
                                X[j,m] = X[j-1,m]
                                Y[j,m] = Y[j-1,m]
                                a[j,m] = a[j-1,m]
                                g[j,m] = g[j-1,m]
                            else:
                                row,col=np.where(dis==np.min(dis))
                                X[j,m] = x[row[0]]
                                Y[j,m] = y[row[0]]
                                a[j,m] = ma[row[0]]
                                g[j,m] = mg[row[0]]
                                
                            if a[j,m]>=10*a[j-1,m] or a[j,m]<=a[j-1,m]/10 or g[j,m]<=g[j-1,m]/10:
                                diss=np.sort(dis,axis=0)
                                if diss[1]>min_d:
                                        X[j,m] = X[j-1,m]
                                        Y[j,m] = Y[j-1,m]
                                        a[j,m] = a[j-1,m]
                                        g[j,m] = g[j-1,m]
                                else:
                                    row, col=np.where(dis==diss[1])
                                    X[j,m] = x[row[0]]
                                    Y[j,m] = y[row[0]]
                                    a[j,m] = ma[row[0]]
                                    g[j,m] = mg[row[0]]
    
            np.save('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/'+str(ui+1)+'cell'+'/X.npy',X)     
            np.save('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/'+str(ui+1)+'cell'+'/Y.npy',Y)
            np.save('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/'+str(ui+1)+'cell'+'/droplet_id.npy',DI)
            
    
    
    ####################################################
    # start tracking multiple cells
    ####################################################        
        X_cluster=np.zeros((np.int32((len(images_name)-1)/2),xyl_cluster.shape[0]))
        Y_cluster=np.zeros((np.int32((len(images_name)-1)/2),xyl_cluster.shape[0]))
        a_cluster=np.zeros((np.int32((len(images_name)-1)/2),xyl_cluster.shape[0]))
        g_cluster=np.zeros((np.int32((len(images_name)-1)/2),xyl_cluster.shape[0]))
        DI_cluster=np.zeros((1,xyl_cluster.shape[0]))
        
        X_cluster[0,:]=xyl_cluster[:,1]
        Y_cluster[0,:]=xyl_cluster[:,0]
        a_cluster[0,:]=xyl_cluster[:,2]
        g_cluster[0,:]=xyl_cluster[:,3]
        DI_cluster[0,:]=xyl_cluster[:,5]
        
        
        
        
        for ui in range(1,uniquexyl.shape[0]-1):                         
            sort_cell=uniquexyl[ui]
            HFCT.make_path('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/'+str(sort_cell)+'cell')
            N1=0
            j=0
            
            for Nt in range (xyl.shape[0]):
                if xyl[Nt,4]==sort_cell:
                    N1=N1+1                
                else:
                    continue
            X=np.zeros((np.int32((len(images_name)-1)/2),N1))
            Y=np.zeros((np.int32((len(images_name)-1)/2),N1))
            Xd=np.zeros((np.int32((len(images_name)-1)/2),N1))
            Yd=np.zeros((np.int32((len(images_name)-1)/2),N1))
            a=np.zeros((np.int32((len(images_name)-1)/2),N1))
            g=np.zeros((np.int32((len(images_name)-1)/2),N1))
            DI=np.zeros((1,N1))
            M=np.zeros((N1,1))
            N2=0
            for Nt in range (xyl.shape[0]):
                if xyl[Nt,4]==sort_cell:
                    N2=N2+1
                    M[N2-1]=Nt
                    X[j,N2-1] = xyl[Nt,0]
                    Y[j,N2-1] = xyl[Nt,1]
                    a[j,N2-1] = xyl[Nt,2]
                    g[j,N2-1] = xyl[Nt,3]
                    DI[j,N2-1] = xyl[Nt,5]
                else:
                    continue
                    
                    
            if N1>0:
                for j in range (1,np.int32((len(images_name)-1)/2)):
                        A=io.imread(os.path.join('C:/Users/data-'+str(field_of_view),images_name[j]))
                        A=A.astype(np.float64)
                        A=A-1200 # subtract certain amount of background noise
                        A[A<0]=0
                        NT=np.uint16(np.std(A))
                        if NT<=300:
                            NT=NT*2
                        bw=extrema.h_maxima(A, NT)
                        L=skimage.measure.label(bw)
                        props=skimage.measure.regionprops(L,A)
                        x=np.zeros((len(props),1))
                        y=np.zeros((len(props),1))
                        ma=np.zeros((len(props),1))
                        mg=np.zeros((len(props),1))
                        for i, properties in enumerate(props):
                            x[i]=properties.centroid[1]
                            y[i]=properties.centroid[0]
                            ma[i]=properties.area
                            mg[i]=np.sum(properties.intensity_image)
                            
                            
                        A_cluster=io.imread(os.path.join('C:/Users/data-'+str(field_of_view),images_name[j]))
                        B_cluster=io.imread(os.path.join('C:/Users/data-'+str(field_of_view),images_name[j+512]))
                        threshold=filters.threshold_otsu(B_cluster)
                        bw_cluster=threshold<B_cluster
                        L_cluster=skimage.measure.label(bw_cluster)
                        props_cluster=skimage.measure.regionprops(L_cluster,A_cluster)
                        x_cluster=np.zeros((len(props_cluster),1))
                        y_cluster=np.zeros((len(props_cluster),1))
                        ma_cluster=np.zeros((len(props_cluster),1))
                        mg_cluster=np.zeros((len(props_cluster),1))
                        for i, properties in enumerate(props_cluster):
                            x_cluster[i]=properties.centroid[1]
                            y_cluster[i]=properties.centroid[0]
                            ma_cluster[i]=properties.area
                            mg_cluster[i]=np.sum(properties.intensity_image)
                        
                        for m in range (xyl_cluster.shape[0]):
                            dis_cluster=np.zeros((len(props_cluster),1))
                            for n in range (len(props_cluster)):
                                dis_cluster[n]=np.sqrt(  np.square(X_cluster[j-1,m]-x_cluster[n]) + np.square(Y_cluster[j-1,m]-y_cluster[n])   )
                            if np.min(dis_cluster)> min_d:
                                X_cluster[j,m] = X_cluster[j-1,m]
                                Y_cluster[j,m] = Y_cluster[j-1,m]
                                a_cluster[j,m] = a_cluster[j-1,m]
                                g_cluster[j,m] = g_cluster[j-1,m]
                            else:
                                row,col=np.where(dis_cluster==np.min(dis_cluster))
                                X_cluster[j,m] = x_cluster[row[0]]
                                Y_cluster[j,m] = y_cluster[row[0]]
                                a_cluster[j,m] = ma_cluster[row[0]]
                                g_cluster[j,m] = mg_cluster[row[0]]
                                
                            if a_cluster[j,m]>=10*a_cluster[j-1,m] or a_cluster[j,m]<=a_cluster[j-1,m]/10 or g_cluster[j,m]<=g_cluster[j-1,m]/10:
                                diss_cluster=np.sort(dis_cluster,axis=0)
                                if diss_cluster[1]>min_d:
                                        X_cluster[j,m] = X_cluster[j-1,m]
                                        Y_cluster[j,m] = Y_cluster[j-1,m]
                                        a_cluster[j,m] = a_cluster[j-1,m]
                                        g_cluster[j,m] = g_cluster[j-1,m]
                                else:
                                    row, col=np.where(dis_cluster==diss_cluster[1])
                                    X_cluster[j,m] = x_cluster[row[0]]
                                    Y_cluster[j,m] = y_cluster[row[0]]
                                    a_cluster[j,m] = ma_cluster[row[0]]
                                    g_cluster[j,m] = mg_cluster[row[0]]
                        
                        for m in range (N1):
                            droplet_label=xyl[M[m].astype(np.int32),5]
                            row,droplet_label_cluster_index=np.where(DI_cluster==droplet_label)
                            cell_cell_dis=np.zeros((xyl.shape[0],1))
                            for k in range(xyl.shape[0]):
                                cell_cell_dis[k]=np.sqrt(  np.square(xyl[M[m].astype(np.int32),0] - xyl[k,0])   +   np.square(xyl[M[m].astype(np.int32),1] - xyl[k,1]  )   )
                            diss_cell_cell_dis=np.sort(cell_cell_dis)
                            if diss_cell_cell_dis[1]>=8:
                                min_d_ns=8
                            else:
                                min_d_ns=np.floor(diss_cell_cell_dis[1])
                                
                            if  droplet_label_cluster_index.shape[0]==1:
                                Xd[j-1,m]=X[j-1,m]+( X_cluster[j,droplet_label_cluster_index]-X_cluster[j-1,droplet_label_cluster_index] )
                                Yd[j-1,m]=Y[j-1,m]+( Y_cluster[j,droplet_label_cluster_index]-Y_cluster[j-1,droplet_label_cluster_index] )
                            elif droplet_label_cluster_index.shape[0]>=2:
                                dis_cell_cluster=np.zeros((droplet_label_cluster_index.shape[0],1))
                                for n in range(droplet_label_cluster_index.shape[0]):
                                    dis_cell_cluster[n]=np.sqrt( np.square(X[0,m]-X_cluster[0,droplet_label_cluster_index[n]]) + np.square(Y[0,m]-Y_cluster[0,droplet_label_cluster_index[n]]) )
                                row_cluster,col_cluster=np.where(dis_cell_cluster==np.min(dis_cell_cluster))
                                Xd[j-1,m]=X[j-1,m]+( X_cluster[j,droplet_label_cluster_index[row_cluster[0]]]-X_cluster[j-1,droplet_label_cluster_index[row_cluster[0]]] )
                                Yd[j-1,m]=Y[j-1,m]+( Y_cluster[j,droplet_label_cluster_index[row_cluster[0]]]-Y_cluster[j-1,droplet_label_cluster_index[row_cluster[0]]] )
                            else:
                                Xd[j-1,m]=X[j-1,m]
                                Yd[j-1,m]=Y[j-1,m]
                            
                            
                            dis=np.zeros((len(props),1))
                            for n in range (len(props)):
                                dis[n]=np.sqrt(  np.square(Xd[j-1,m]-x[n]) + np.square(Yd[j-1,m]-y[n])   )
                            if np.min(dis)> min_d_ns:
                                X[j,m] = X[j-1,m]
                                Y[j,m] = Y[j-1,m]
                                a[j,m] = a[j-1,m]
                                g[j,m] = g[j-1,m]
                            else:
                                row,col=np.where(dis==np.min(dis))
                                X[j,m] = x[row[0]]
                                Y[j,m] = y[row[0]]
                                a[j,m] = ma[row[0]]
                                g[j,m] = mg[row[0]]
                                
                            if a[j,m]>=10*a[j-1,m] or a[j,m]<=a[j-1,m]/10 or g[j,m]<=g[j-1,m]/10:
                                diss=np.sort(dis,axis=0)
                                if diss[1]>min_d_ns+2:
                                        X[j,m] = X[j-1,m]
                                        Y[j,m] = Y[j-1,m]
                                        a[j,m] = a[j-1,m]
                                        g[j,m] = g[j-1,m]
                                else:
                                    row, col=np.where(dis==diss[1])
                                    X[j,m] = x[row[0]]
                                    Y[j,m] = y[row[0]]
                                    a[j,m] = ma[row[0]]
                                    g[j,m] = mg[row[0]]
    
            np.save('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/'+str(sort_cell)+'cell'+'/X.npy',X)     
            np.save('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/'+str(sort_cell)+'cell'+'/Y.npy',Y)
            np.save('C:/Users/data-results/local maxium separate cell/'+str(field_of_view)+'/'+str(sort_cell)+'cell'+'/droplet_id.npy',DI)
    return 0                                        
                        