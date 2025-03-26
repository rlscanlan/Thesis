#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 15:34:16 2022

@author: rebekahscanlan
"""

from pathlib import Path
import os, glob
import pandas, numpy
from collections import OrderedDict
import tellurium as te
import site
from sys import platform
import matplotlib.pyplot as plt
import seaborn
import matplotlib
from _simulator import TimeSeriesPlotter, TimeSeries

# from data.data_analysis import
from TGFmodel import *
import site

import teUtils

import re 

matplotlib.use('TkAgg')
seaborn.set_context('talk')


class TheModel:
    model_string="""
    function MM(km, Vmax, S)
        Vmax * S / (km + S)
    end
    
    function MMWithKcat(km, kcat, S, E)
        kcat * E * S / (km + S)
    end
    
    function NonCompetitiveInhibition(km, ki, Vmax, I, S)
        Vmax * S / ((km + S) * (1 + I / ki))
    end
    
    function NonCompetitiveInhibitionWithKcat(km, ki, kcat, E, I, S)
        kcat * E * S / ( (km + S) * (1 + (I / ki) ) )
    end
    
    function NonCompetitiveInhibitionWithKcatAndExtraActivator(km, ki, kcat, E1, E2, n, I, S)
        kcat * E1 * E2 * S / ( (km + S) * (1 + (I / ki)^n ) )
    end
    
    function NonCompetitiveInhibitionWithKcatAndExtraActivatorNoSub(km, ki, kcat, E1, E2, n, I)
        kcat * E1 * E2 / ( km * (1 + (I / ki)^n ) )
    end
    
    function CompetitiveInhibitionWithKcat(km, ki, kcat, E, I, S)
        kcat * E * S / (km + S + ((km * I )/ ki)  )
    end
    
    function CompetitiveInhibition(Vmax, km, ki, I, S)
        Vmax * S / (km + S + ((km * I )/ ki)  )
    end

    function NonCompetitiveInhibitionNoSub(Vmax, km, ki, I)
        Vmax / (km * (1 + I /ki))
    end

    function NonCompetitiveInhibitionWithKcatNoSub(km, ki, kcat, E, I)
        kcat * E / (km * (1 + (I / ki) ))
    end

    function Hill(km, kcat, L, S, h)
        kcat * L * (S / km)^h  /   1 + (S / km)^h
    end
    
    model SimpleTGF_SenModel()
        compartment                 Cell = 1;
        //Cell cycle
        var p53                     in Cell;
        var pp53                    in Cell;
        var ppRb                    in Cell;
        var pRb                     in Cell;
        var p21                     in Cell;
        var p16                     in Cell;
        var DDR                     in Cell;
        //p38
        var p38                     in Cell;
        var pp38                    in Cell;
        //Ras
        var RAS                     in Cell;
        
        // insert inputs here

        DDR = 5;
        RAS = 5;

        // insert initial concentration parameters here:
        //Cell cycle arrest
        p53                     = 1;
        pp53                    = 0.1;
        ppRb                    = 1.25;
        pRb                     = 0.05;
        p21                     = 0.1;
        p16                     = 0.1;
        DDR                     = 5;

        //p38
        p38                     = 1;
        pp38                    = 0.1;

        //Metabolism
        RAS                     = 5;
        
        //Kinetic parameters
        
        //Cell Cycle Arrest
        kp53A                           = 0.5;
        kpp53B                          = 0.5;
        kp53D                           = 5;
        kp53F                           = 0.75;
        kppRbA                          = 0.25;
        VmaxppRbB                       = 0.25;
        kppRbB                          = 0.2;
        kInppRbBByp16                   = 0.3;
        kInppRbBByp21                   = 0.3;
        kp21F                           = 1;
        kp21D                           = 0.5;
        kp16F1                          = 0.3;
        kp16F2                          = 0.5;
        kp16D                           = 0.25;
        kDDRD                           = 6;
        kDDRFRAS                        = 0.5;
        kDDRF                           = 0.5;

        //p38
        kp38A                           = 0.4;    
        kpp38B                          = 0.3;
        kp38F                           = 0.7;
        kp38D                           = 0.5;

        //Ras
        kRASD                           = 0.2;
        kRASF                           = 0.1;

        n                               = 1;   
        
           //Cell cycle arrest
       CR1F3    :   => p21                                ; Cell * kp21F * pp53;
       CR1D1    :   p21 =>                                ; Cell * kp21D * p21;
       CR3F1    :   => p53                                ; Cell * kp53F * DDR;
       CR3D1    :   p53 =>                                ; Cell * kp53D * p53;
       CR3A1    :   p53 => pp53                           ; Cell * kp53A * p53 * DDR;
       CR3B1    :   pp53 => p53                           ; Cell * kpp53B * pp53;
       CR4F1    :   => p16                                ; Cell * kp16F1 * pp38;
       CR4F2    :   => p16                                ; Cell * kp16F2 * DDR;
       CR4D1    :   p16 =>                                ; Cell * kp16D * p16;
       CR5B1    :   pRb => ppRb                           ; Cell * NonCompetitiveInhibition(kppRbB, kInppRbBByp21, VmaxppRbB, p21, pRb);
       CR5B2    :   pRb => ppRb                           ; Cell * NonCompetitiveInhibition(kppRbB, kInppRbBByp16, VmaxppRbB, p16, pRb);
       CR5A1    :   ppRb => pRb                           ; Cell * kppRbA * ppRb;
       CR6F1    :   => DDR                                ; Cell * kDDRF;
       CR6F2    :   => DDR                                ; Cell * kDDRFRAS * RAS;
       CR6D1    :   DDR =>                                ; Cell * kDDRD * DDR;
        //p38 dynamics
       IR6B1    :   pp38 => p38                           ; Cell * kpp38B * pp38;
       IR6A1    :   p38 => pp38                           ; Cell * kp38A * p38 * RAS;
       IR6F1    :   => p38                                ; Cell * kp38F;
       IR6D1    :   p38 =>                                ; Cell * kp38D * p38;
        //RAS 
       MR3F1    :   => RAS                                ; Cell * kRASF;
       MR3D1    :   RAS =>                                ; Cell * kRASD * RAS * DDR;
     
       //Events
       //DDIS
       at 35 after (time>0): DDR = 5;
       at 35 after (time>0): kDDRF = 14;
       at 35 after (time>0): kRASF = 0.3;
       //OIS
       //at 35 after (time>0): RAS = 5;
       //at 35 after (time>0): kDDRFRAS = 5.75;
       //at 35 after (time>0): kRASF = 1.75;
       //p53 knockdown T35
       at 35 after (time>0): p53 = 0.01;
       at 35 after (time>0): kp53F = 0.1;
       //at 35 after (time>0): kp53F2 = 0.1;

    end
    
    """

    def __init__(self):
        self.rr = self._load_rr()

    def _load_rr(self):
        return te.loada(self.model_string)


    def simulate(self, start, stop, num):
        data = self.rr.simulate(start, stop, num)
        data = pandas.DataFrame(data, columns=data.colnames).set_index('time')
        data.columns = [i.replace('[', '').replace(']', '') for i in data.columns]
        return data

if __name__ == '__main__':

    OPEN_WITH_COPASI = False

    PLOT_SIMULATION = True

    PARAMETER_SCAN = True


    # ========================================

    mod = TheModel()

    if PLOT_SIMULATION:
    
        NonSen = {
            'DDR': 0,
            'RAS': 0.5,
            'kDDRF': 0.01,
            'kRASF': 0.01,
            #'kDDRFRAS': 1,
        }
    
        plot_selection = {
            'Model A - Ras and p38': ['RAS'],
            'Model A - Cell cycle_1': ['p21', 'p16','DDR'],
            'Model A - Cell cycle_2': ['p53', 'pp53', 'ppRb', 'pRb'],
            'Model A - p38': ['p38', 'pp38'],
        }

        ts = TimeSeries(mod.model_string, NonSen, 0, 100, 1000)
        tsp = TimeSeriesPlotter(ts, plot_selection, plot_dir=SIMULATION_DIRECTORY, fname='DDIS_A_p53KD.png',
                          savefig=True, legend=True, ncols=3, ylim=(-0.5,10))
        tsp.plot()
        tsp.data.to_csv('/Users/rebekahscanlan/python/tellurium/TGFmodel/data_files/DDIS_A_p53KD.csv')
        print(ts)

