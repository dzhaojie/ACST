# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 22:22:57 2017

@author: dzhaojie
"""

import ExtractIntensityMulti as EIM
import multiprocessing
import numpy as np




if __name__ == '__main__':
    

################################
#process all the field of views
################################

   
   jobs = []
   for i_num_of_field in range(12):
        num_of_field=range(i_num_of_field*7,i_num_of_field*7+7)
        p = multiprocessing.Process(target=EIM.extract_intensity, args=(num_of_field,))
        jobs.append(p)
        p.start()
    
   for i_j in jobs:
        i_j.join()
       

################################
#process certain field of views
################################

   jobs = []
   fov=list(np.array([82]))
   for i_num_of_field in range(1):
        num_of_field=fov[i_num_of_field*1:i_num_of_field*1+1]
        p = multiprocessing.Process(target=EIM.extract_intensity, args=(num_of_field,))
        jobs.append(p)
        p.start()
    
   for i_j in jobs:
        i_j.join()

        