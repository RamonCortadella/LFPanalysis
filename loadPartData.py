#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 17:19:03 2022

@author: ramon
"""
import h5py
import numpy as np
import gc
import quantities as pq
import deepdish as dd
from PhyREC.NeoInterface import NeoSegment, NeoSignal
path = '../../../data/LargeScale/'
FilesList = [
            
            'B13289O14-DH1-0164/Day1-09_10-12-21/RawData/B13289O14-DH1-Rec1-v2good-Mocap-Vgs0p2V-Vds0p1V_0.h5',
            'B13289O14-DH1-0164/Day1-09_10-12-21/RawData/B13289O14-DH1-Rec2-Mocap-Vgs0p17V-Vds0p1V_0.h5',
            'B13289O14-DH1-0164/Day1-09_10-12-21/RawData/B13289O14-DH1-Rec3-Mocap-Vgs0p13V-Vds0p1V_0.h5',
            'B13289O14-DH1-0164/Day1-09_10-12-21/RawData/B13289O14-DH1-Rec4-v2good-Mocap-Vgs0p11V-Vds0p1V_0.h5',
            'B13289O14-DH1-0164/Day1-09_10-12-21/RawData/B13289O14-DH1-Rec5-Mocap-Vgs0p11V-Vds0p1V_0.h5',
            'B13289O14-DH1-0164/Day1-09_10-12-21/RawData/B13289O14-DH1-Rec6-v2-Mocap-minusVgs0p05V-Vds0p1V_0.h5',
            'B13289O14-DH1-0164/Day1-09_10-12-21/RawData/B13289O14-DH1-Rec7-v2good-Mocap-Vgs0p0V-Vds0p1V_0.h5',
            'B13289O14-DH1-0164/Day1-09_10-12-21/RawData/B13289O14-DH1-Rec8-Mocap-Vgs0p05V-Vds0p1V_0.h5',
            'B13289O14-DH1-0164/Day1-09_10-12-21/RawData/B13289O14-DH1-Rec9-v2ggood-Mocap-minVgs0p025V-Vds0p1V_0.h5',
            
            
            'B13289O14-DH2-01551/RawData/B13289O14-DH2-Rec1-v3-Vgs0p28-Vds0p1_0.h5',
            'B13289O14-DH2-01551/RawData/B13289O14-DH2-Rec2-Vgs0p28-Vds0p1_0.h5',
            'B13289O14-DH2-01551/RawData/B13289O14-DH2-Rec3-Vgs0p28-Vds0p1_0.h5',
            'B13289O14-DH2-01551/RawData/B13289O14-DH2-Rec4-Vgs0p28-Vds0p1_0.h5',
            'B13289O14-DH2-01551/RawData/B13289O14-DH2-Rec11-Vgs0p28-Vds0p1_0.h5',
            'B13289O14-DH2-01551/RawData/B13289O14-DH2-Rec12-Vgs0p28-Vds0p1_0.h5',
            'B13289O14-DH2-01551/RawData/B13289O14O23-DH2SL4-Rec1-Vgs0p25-Vds0p1_0.h5',
            
            'B13289O14-DH3-01553/RawData/B13289O14-DH3-Rec1-Vgs0p3Vds0p1_0.h5',
            'B13289O14-DH3-01553/RawData/B13289O14-DH3-Rec2-Vgs0p3Vds0p1_0.h5',
            'B13289O14-DH3-01553/RawData/B13289O14-DH3-Rec3-v2-Vgs0p3Vds0p1_0.h5',
            'B13289O14-DH3-01553/RawData/B13289O14O23-DH3SL5-Rec1-v2-Vgs0p3Vds0p1_0.h5',
            'B13289O14-DH3-01553/RawData/B13289O14O23-DH3SL5-Rec2-Vgs0p32Vds0p1_0.h5',
            'B13289O14-DH3-01553/RawData/B13289O14O23-DH3SL5-Rec3-Vgs0p32Vds0p1_0.h5',
            'B13289O14-DH3-01553/RawData/B13289O14O23-DH3SL5-Rec4-v4-Vgs0p24Vds0p1_0.h5',
            'B13289O14-DH3-01553/RawData/B13289O14O23-DH3SL5-Rec5-Vgs0p26Vds0p1_0.h5',
            'B13289O14-DH3-01553/RawData/B13289O14O23-DH3SL5-Rec6-Vgs0p22Vds0p1_0.h5',
            
            'B13289O14-DH5-01603/RawData/B13289O14-DH5-Rec1-Vgs0p36-Vds0p1_0.h5',
            'B13289O14-DH5-01603/RawData/B13289O14-DH5-Rec2-Vgs0p36-Vds0p1_0.h5',
            'B13289O14-DH5-01603/RawData/B13289O14-DH5-Rec3-v2-Vgs0p36-Vds0p1_0.h5',
            'B13289O14-DH5-01603/RawData/B13289O14-DH5-Rec4-Vgs0p36-Vds0p1_0.h5',
            
            'B13289O24-DH1-01604/RawData/B13289O24-DH1SL7-Rec3-v5-Vgsm0p35-ElectFixed-Vds0p1_0_0.h5',
            'B13289O24-DH1-01604/RawData/B13289O24-DH1SL7-Rec4-Vgsm0p35-ElectFixed-Vds0p1_0_0.h5',
            'B13289O24-DH1-01604/RawData/B13289O24-DH1SL7-Rec5-Vgsm0p3-ElectFixed-Vds0p1_0_0.h5',
            
            ]
for File in FilesList:
    with h5py.File(path+File,'r') as h5f:
        DataRows = h5f['data'][:,1]
        
    gc.collect()
    nr = 32
    ncols = 16
    
    Fs = 1e6*pq.Hz
    SampsSwitchPeriod = 3
    StabTimeSw = 20e-6 #in 
    
    FsCh = Fs/(nr*SampsSwitchPeriod*ncols)
    
    SSP = SampsSwitchPeriod
    nc = ncols
    
    CroppedDR = DataRows.reshape((SSP, -1), order='F')
    
    ic = 0;
    
    chdat = CroppedDR[:,ic::nc]
    StabInd = int(np.floor(StabTimeSw*Fs/(nr*pq.Hz)))+1
    s = np.mean(chdat[-StabInd:,:],axis=0)
    
    
    sig = NeoSignal(s,
                    units='A',
                    sampling_rate=FsCh,
                    name='Ch01Col01DC')
    f = File.split('/')
    fN = f[0].split('-')
    path2 = path + f[:-2]
    FilePathAndName = path2 + 'CalibrationScripts/'+ '-'.join(fN[1:4])+'/'+'triggerSignal.h5'
    print(FilePathAndName)
    dd.io.save(FilePathAndName, sig)