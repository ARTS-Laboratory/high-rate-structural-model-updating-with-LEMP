# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 17:35:05 2024

@author: Alexander
"""
import os
import numpy as np
from collections import namedtuple
from collections import defaultdict

sm = namedtuple('Signal_Metrics', ['signal','TRAC','SNR','meanErr', 'errVar', 'errMax','errMin','samples'])
tm = namedtuple('Timinig_Metrics', ['timed','mean', 'var', 'max','min','samples'])
ts = namedtuple('Time_Series', ['Time','Meas_Pos','Model_Pos', 'Filter_Pos','Meas_Freq','Model_Freq', 'Filter_Freq'])
labels = ["FEANodes","Dataset","ResampleMethod","EigenSolver","ModelsForcasted","FreqKF","StateKF","BufferSizeTimeSeries"]
tags = ["Node","Data","Samp","Eign","Mods","FrKF","StKF","Buff"]
trial_props = namedtuple('Trail_Details', tags)




def list_files_in_directory(directory):
    stats_files = []
    time_series_files = []
    for filename in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, filename)) and 'statistics'in filename:
            stats_files.append(filename)
        if os.path.isfile(os.path.join(directory, filename)) and 'time_series'in filename:
            time_series_files.append(filename)
    return stats_files, time_series_files




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
    def __init__(self,trail_details,time_series:str,
                 signal_metrics:dict={},timinig_metrics:dict={})->None:
        self.signal_metrics=signal_metrics
        self.timinig_metrics=timinig_metrics
        self.trail_details=trail_details
        self.time_series=time_series
        self.time_series_data = None
        
    def setter(self,trail_details,time_series:str,
               signal_metrics:dict={},timinig_metrics:dict={})->None:
        self.signal_metrics=signal_metrics
        self.timinig_metrics=timinig_metrics
        self.trail_details=trail_details
        self.time_series=time_series
        
    def genKey(self)->str:
        p=self.trail_details
        return f"{p.Node},{p.Samp},{p.Eign},{p.Mods}"
    
    def loadTimeSeries(self)->tuple:
        if self.time_series_data == None:
            data = np.loadtxt(self.time_series,skiprows=23,delimiter=',')
            self.time_series_data = ts(*[data[:,i] for i in range(np.shape(data)[1])])
        else: 
            pass
        return self.time_series_data

def loadExperiments():
    current_directory = os.getcwd()
    files = list_files_in_directory(current_directory)
    experiments = defaultdict(lambda : False)
    
    for sfile,tsfile in zip(files[0],files[1]):
        with open(sfile, 'r') as stat:
            ex=experiment(trial_props, tsfile)
            for i, line in enumerate(stat):
                line_sl=line.split(",")
                if "Signal Characteristics" in line_sl[0]:
                    data=[float(d) for d in line_sl[1:]]
                    data.insert(0,line_sl[0])
                    ex.signal_metrics.update({line_sl[0]:sm(*data)})
                    
                if "Timing Distribution" in line_sl[0]:
                    data=[float(d) for d in line_sl[1:]]
                    data.insert(0,line_sl[0])
                    ex.timinig_metrics.update({line_sl[0]:tm(*data)})
    
                if line_sl[0]== labels[0]:
                    ex.trail_details.Node=line_sl[1].strip('\n')
                elif line_sl[0]== labels[1]:
                    ex.trail_details.Data=line_sl[1].strip('\n')
                elif line_sl[0]== labels[2]:
                    ex.trail_details.Samp=line_sl[1].strip('\n')
                elif line_sl[0]== labels[3]:
                    ex.trail_details.Eign=line_sl[1].strip('\n')
                elif line_sl[0]== labels[4]:
                    ex.trail_details.Mods=line_sl[1].strip('\n')
                elif line_sl[0]== labels[5]:
                    ex.trail_details.FrKF=line_sl[1].strip('\n')
                elif line_sl[0]== labels[6]:
                    ex.trail_details.StKF=line_sl[1].strip('\n')
                elif line_sl[0]== labels[6]:
                    ex.trail_details.Buff=line_sl[1].strip('\n')
                else:
                    pass
        experiments.update({ex.genKey():experiment})
    return experiments

Experiments=loadExperiments()
def fx(x):
    return x**2 +2*x +1
a=2
function_example = np.vectorize(lambda x: x**2 +a*x +1)
x=np.linspace(0,100,num=1000)
y=function_example(x)


