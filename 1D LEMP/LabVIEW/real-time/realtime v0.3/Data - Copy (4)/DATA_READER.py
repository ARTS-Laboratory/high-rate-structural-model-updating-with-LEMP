# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 16:37:21 2024

@author: Alexander
"""
import numpy as np
import os


gaussian = lambda x, u, s: 1/(np.sqrt(2*np.pi)*sig)*np.exp(-np.power((x-mu)/sig,2)/2)
CI95 = lambda mu, s, n: mu+1.96*(s/(n**.5)), mu-1.96*(s/(n**.5)) 
normal_dist = np.vectorize(gaussian)


class guassian_stat:
    def __init__(self, mean:float=0,variance:float=0,maximum:float=0,
                 minimum:float=0,samples:int=0):
        set_data(self, mean,variance,maximum,minimum,samples)
        return self
    
    def set_data(self, mean:float=0,variance:float=0,maximum:float=0,
                 minimum:float=0,samples:int=0):
        self.mean,self.variance,self.maximum,self.minimum,self.samples=mean,variance,maximum,minimum,samples
        return self
    
    def generate_distribution(self, length:int=400, bounds:iter=[0,4000]):
        mu,d_max,d_min=self.mean,self.maximum, self.minimum
        s,n=np.sqrt(self.variance),self.samples
        X=np.linspace(bounds[0], bounds[1], num=length)
        gaussian = np.vectorize(
            lambda x : 1/(np.sqrt(2*np.pi)*s)*np.exp(-np.power((x-mu)/s,2)/2))
        dist = gaussian(X)
        dist = (dist-np.min(dist))/(np.max(dist)-np.min(dist))
        ci=CI95(mu,sig,n)
        plot_min_max =([1,0],[d_max,d_max],[d_min,d_min])
        
        return 
        
        
        
    
class signal_stat:
    def __init__(self, TRAC:float=0, SNR:float=0, mean_err:float=0,
                 variance_err:float=0,maximum_err:float=0,minimum_err:float=0,
                 samples:int=0):
        set_data(self, mean,variance,maximum,minimum,samples)
        return self
    
    def set_data(self, TRAC:float=0, SNR:float=0, mean_err:float=0,
                 variance_err:float=0,maximum_err:float=0,minimum_err:float=0,
                 samples:int=0):
        self.mean,self.variance,self.maximum,self.minimum,self.samples=mean,variance,maximum,minimum,samples
        return self

class experiment_statistics:
    def __init__(self, properties:dict={},timing_statistics:dict={},
                 signal_characteristics:dict={}):
        set_data(self, properties,timing_statistics,signal_characteristics)
        return self
    
    def set_data(self, properties:dict={},timing_statistics:dict={},
                 signal_characteristics:dict={}):
        self.properties,self.timing_statistics,self.signal_characteristics=properties,timing_statistics,signal_characteristics
        return self
    
    def add_timing_distribution(self, name:str='',mean:float=0,
                                variance:float=0,maximum:float=0,
                                minimum:float=0,samples:int=0):
        self.timing_statistics.update({name:guassian_stat(mean,variance,maximum,minimum,samples)})
        return self
    
    def add_properties(self, properties:dict={}):
        self.properties.update(properties)
        return self
    
    def add_properties(self, properties:dict={}):
        self.properties.update(properties)
        return self
    
    

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

current_directory = os.getcwd()
files = list_files_in_directory(current_directory)
print("Files in directory:")
for file in files:
    print(file)
    with open(file):
        
    
