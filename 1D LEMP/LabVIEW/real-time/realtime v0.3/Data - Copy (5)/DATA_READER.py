# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 16:37:21 2024

@author: Alexander
"""
import numpy as np
import os



#Norm_Series = lambda X , F :



class guassian_stat:
    def __init__(self,mean:float=0,variance:float=0,maximum:float=0,
                 minimum:float=0,samples:int=0):
        self.set_data(self, mean,variance,maximum,minimum,samples)
        return self
    
    def set_data(self,mean:float=0,variance:float=0,maximum:float=0,
                 minimum:float=0,samples:int=0):
        self.mean,self.variance,self.std_dev=mean,variance,np.sqrt(variance)
        self.maximum,self.minimum,self.samples=maximum,minimum,samples
        return self
    
    def gen_distribution(self,length:int=400,bounds:iter=[0,4E3],norm:bool=True):
        mu,d_max,d_min=self.mean,self.maximum, self.minimum
        s,n=np.sqrt(self.variance),self.samples
        CI95 = lambda s, n: 1.96*(s/(n**.5))
        gaussian = np.vectorize(
            lambda x : 1/(np.sqrt(2*np.pi)*s)*np.exp(-np.power((x-mu)/s,2)/2))
        X=np.linspace(bounds[0], bounds[1], num=length)
        
        dist=gaussian(X)
        if norm:
            dist=(dist-np.min(dist))/(np.max(dist)-np.min(dist))
        ci0, ci1= mu+CI95(s,n),mu-CI95(s,n)
        plot_thresh=([1,0],[d_max,d_max],[d_min,d_min],[ci0,ci0],[ci1,ci1])        
        return dist, plot_thresh
    
    def get_plottables(self,length:int=400,bounds:iter=[0,4E3],norm:bool=True):
        return self.gen_distribution(length,bounds,norm)

        
        
    
class signal_stat(guassian_stat):
    def __init__(self, TRAC:float=0,SNR:float=0,mean:float=0,
                 variance:float=0,maximum:float=0,minimum:float=0,
                 samples:int=0):
        self.set_data(self,TRAC,SNR,mean,variance,maximum,minimum,samples)
        return self
    
    def set_data(self, TRAC:float=0,SNR:float=0,mean:float=0,
                 variance:float=0,maximum:float=0,minimum:float=0,
                 samples:int=0):
        self.TRAC,self.SNR=TRAC,SNR
        self.mean_err,self.variance_err,self.std_dev_err=mean,variance,np.sqrt(variance)
        self.maximum_err,self.minimum_err,self.samples=maximum,minimum,samples
        return self

    def get_plottables(self,length:int=400,bounds:iter=[0,4E3],norm:bool=True):
        return self.gen_distribution(length,bounds,norm),[self.TRAC,self.SNR]
    

    
class experiment:
    def __init__(self, stats_file:str, time_series_file:str):
        self.sFname, self.tsFname = stats_file, time_series_file
        
    def load_stats(self):
        with open(self.sfname) as stats:
            fSignal = stats[2].split(",")[1:]
            sSignal = stats[3].split(",")[1:]
            pass
    
    def load_time_series(self):

# class experiment_statistics:
#     def __init__(self, properties:dict={},timing_statistics:dict={},
#                  signal_characteristics:dict={}):
#         self.set_data(self, properties,timing_statistics,signal_characteristics)
#         return self
    
#     def set_data(self, properties:dict={},timing_statistics:dict={},
#                  signal_characteristics:dict={}):
#         self.properties,self.timing_statistics,self.signal_characteristics=properties,timing_statistics,signal_characteristics
#         return self
    
#     def add_timing_distribution(self, name:str='',mean:float=0,
#                                 variance:float=0,maximum:float=0,
#                                 minimum:float=0,samples:int=0):
#         self.timing_statistics.update({name:guassian_stat(mean,variance,maximum,minimum,samples)})
#         return self
    
#     def add_properties(self, properties:dict={}):
#         self.properties.update(properties)
#         return self
    

    

# def file_next(base_name, extension, increment=1, start=1, flist=set()):
#     index = start
#     while True:
#         file_name = f"{base_name}{index}.{extension}"
#         if os.path.exists(file_name):
#             flist.add(file_name)
#         else:
#             return flist
#         index += increment

def list_files_in_directory(directory):
    files = []
    for filename in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, filename)):
            files.append(filename)
    return files

# current_directory = os.getcwd()
# files = list_files_in_directory(current_directory)
# print("Files in directory:")
# for file in files:
#     print(file)
#     with open(file):
        
    
