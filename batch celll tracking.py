# -*- coding: utf-8 -*-
"""
Created on Sat Oct 21 23:24:45 2017

@author: dzhaojie
"""
import CellTrackingMulti as CTM
import multiprocessing




if __name__ == '__main__':
    
#    num_of_fields=range(84)
#    #num_cores = multiprocessing.cpu_count()
#    pool = multiprocessing.Pool(processes=3)
#    [pool.apply(CTM.cell_track, args=(num_of_field,)) for num_of_field in num_of_fields]
#    pool.start()
    
   jobs = []
   for i_num_of_field in range(12):
        num_of_field=range(i_num_of_field*7,i_num_of_field*7+7)
        p = multiprocessing.Process(target=CTM.cell_track, args=(num_of_field,))
        jobs.append(p)
        p.start()
    
   for i_j in jobs:
        i_j.join()
      
    