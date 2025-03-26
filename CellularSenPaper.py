#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 09:55:47 2022

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
from libsbml import *
import simplesbml


# from data.data_analysis import
from TGFmodel import *
import site

import teUtils

from pycotools3 import tasks, viz, model

import re 

def addCopasiPath(copasiPath):
    """adds the path to copasi to pythons working path

    Copasi is often not in the path python uses. This will check if its
    not and add it.

    Args:
        copasiPath (str): a string with the path to copasi on the
        relivent machine.
    """
    if not re.search(copasiPath, os.environ["PATH"]):
       os.environ["PATH"] += os.pathsep + copasiPath
       
addCopasiPath("/Applications/COPASI")   

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
        //pInflamm & p38
        var p38                     in Cell;
        var pp38                    in Cell;
        var p65IkBa                 in Cell;
        var NFkB                    in Cell;
        var pNFkB                   in Cell;
        var IL6                     in Cell;
        var IL1a                    in Cell;
        var CEBPb                   in Cell;
        //Ras
        var RAS                     in Cell;
        var PTEN                    in Cell;
        var mTOR                    in Cell;
        //Notch
        var NOTCH1                  in Cell;
        var NICD                    in Cell;

        // insert initial concentration parameters here:
        //Cell cycle arrest
        p53                     = 1;
        pp53                    = 0.1;
        ppRb                    = 1.25;
        pRb                     = 0.05;
        p21                     = 0.1;
        p16                     = 0.1;
        DDR                     = 5;

        //Inflammation and p38
        IL1a                    = 0.1;
        CEBPb                   = 0.1;
        p65IkBa                 = 1;
        NFkB                    = 0.1;
        pNFkB                   = 0.1;
        IL6                     = 0.1;
        p38                     = 1;
        pp38                    = 0.1;

        //Metabolism
        RAS                     = 5;
        PTEN                    = 1;
        mTOR                    = 1;
        
        //Notch signalling
        NICD                    = 0.1;
        NOTCH1                  =1;
        
        //Kinetic parameters
        
        //Cell Cycle Arrest
        kp53A                           = 0.5;
        kpp53B                          = 0.5;
        kp53D                           = 5;
        kp53F                           = 0.75;
        kp53F2                          = 0.15;
        Vmaxp53F                        = 0.4;
        kp53F3                          = 0.2;
        kp53IbyRAS                      = 0.1;
        kppRbA                          = 0.25;
        VmaxppRbB                       = 0.25;
        kppRbB                          = 0.2;
        kInppRbBByp16                   = 0.3;
        kInppRbBByp21                   = 0.3;
        kp21F1                          = 0.5;
        kp21F2                          = 1;
        kp21D                           = 0.5;
        kp16F1                          = 0.3;
        kp16F2                          = 0.3;
        kp16D                           = 0.25;
        kDDRD                           = 6;
        kDDRFRAS                        = 0.5;
        kDDRF                           = 0.5;

        //Inflammation and p38
        kp65A                           = 0.5;
        kNFkBB                          = 0.4;
        kNFkBA                          = 0.6;
        kpNFkBB                         = 0.4;
        kIL6FbypNFkB                    = 0.4;
        kIL6FIbypp53                    = 1.5;                
        kcatIL6F                        = 0.3;
        kIL6F                           = 0.8;
        kIL6D                           = 0.8;
        kp65F                           = 0.5;
        kp65D                           = 0.5;
        kp38A                           = 0.1;  
        kpp38B                          = 0.25; 
        kp38ARas                        = 0.25;
        kp38AInByNICD                   = 0.1;
        kcatp38ARas                     = 0.1;
        kp38F                           = 0.7;
        kp38D                           = 0.5;
        kIL1aF1                         = 0.5;
        kIL1aF2                         = 0.2;
        kIL1aD                          = 0.5;
        VmaxCEBPbF                      = 0.5;
        kCEBPbF                         = 0.75;
        kCEBPbIn                        = 0.125;
        kCEBPbD                         = 0.5;  

        //Ras
        kRASD                           = 0.2;
        kRASF                           = 0.1;
        kPTENF                          = 0.5;
        kPTEND                          = 0.5;
        kmTORF                          = 0.5;
        kmTORD                          = 5;

        
        //Notch signalling
        kNOTCH1F1                       = 1;
        kNOTCH1F2                       = 0.125;
        kNOTCH1D                        = 0.5;
        kNOTCH1A                        = 1;
        kNICDD                          = 0.75;

        n                               = 1;   
        
        //Cell cycle arrest
       CR1F1    :   => p21                                ; Cell * kp21F1 * NICD;
       CR1F3    :   => p21                                ; Cell * kp21F2 * pp53;
       CR1D1    :   p21 =>                                ; Cell * kp21D * p21;
       CR3F1    :   => p53                                ; Cell * kp53F * DDR;
       CR3F2    :   => p53                                ; Cell * kp53F2 * pNFkB;
       CR3F3    :   => p53                                ; Cell * NonCompetitiveInhibitionNoSub(Vmaxp53F, kp53F3, kp53IbyRAS, RAS);  
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
        //Inflammation and p38 dynamics
       IR1F1   :   => p65IkBa                            ; Cell * kp65F;
       IR1D1   :   p65IkBa =>                            ; Cell * kp65D * p65IkBa;
       IR1A1   :   p65IkBa => NFkB                       ; Cell * kp65A * p65IkBa * IL1a;
       IR1B1   :   NFkB => p65IkBa                       ; Cell * kNFkBB * NFkB;
       IR2F1   :   NFkB => pNFkB                         ; Cell * kNFkBA * NFkB * pp38; 
       IR2B1   :   pNFkB => NFkB                         ; Cell * kpNFkBB * pNFkB;
       IR4D1   :   IL6 =>                                ; Cell * kIL6D * IL6; 
       IR4F1   :   => IL6                                ; Cell * kIL6F * pNFkB;
       IR6B1   :   pp38 => p38                           ; Cell * kpp38B * pp38;
       IR6A1   :   p38 => pp38                           ; Cell * kp38A * p38;
       IR6A2   :   p38 => pp38                           ; Cell * NonCompetitiveInhibitionWithKcatNoSub(kp38ARas, kp38AInByNICD, kcatp38ARas, RAS, NICD); 
       IR6F1   :   => p38                                ; Cell * kp38F;
       IR6D1   :   p38 =>                                ; Cell * kp38D * p38; 
       IR3F1   :   => IL1a                               ; Cell * kIL1aF1 * CEBPb;
       IR3F2   :   => IL1a                               ; Cell * kIL1aF2;  
       IR3D1   :   IL1a =>                               ; Cell * kIL1aD * IL1a;
       IR5F1   :   => CEBPb                              ; Cell * NonCompetitiveInhibitionNoSub(VmaxCEBPbF, kCEBPbF, kCEBPbIn, NICD);
       IR5D1   :   CEBPb =>                              ; Cell * kCEBPbD * CEBPb;       
        //RAS
       MR3F1   :   => RAS                                ; Cell * kRASF;
       MR3D1   :   RAS =>                                ; Cell * kRASD * RAS * DDR;    
       //Notch signalling
       NR1F1   :   => NOTCH1                             ; Cell * kNOTCH1F1 * pp53;
       NR1F2   :   => NOTCH1                             ; Cell * kNOTCH1F2;
       NR1D1   :   NOTCH1 =>                             ; Cell * kNOTCH1D * NOTCH1;
       NR1A1   :   NOTCH1 => NICD                        ; Cell * kNOTCH1A * NOTCH1;
       NR2D1   :   NICD =>                               ; Cell * kNICDD * NICD;
     
       //Events
       //DDIS
       at 35 after (time>0): DDR = 5;
       at 35 after (time>0): kDDRF = 14;
       at 35 after (time>0): kRASF = 0.3;
       //OIS
       //at 35 after (time>0): RAS = 5;
       //at 35 after (time>0): kDDRFRAS = 5.75;
       //at 35 after (time>0): kRASF = 1;
       //Notch switch
       at 65 after (time>0): NICD = 0.01;
       at 65 after (time>0): kNICDD = 5;
       //p53 knockdown T35
       //at 35 after (time>0): p53 = 0.01;
       //at 35 after (time>0): kp53F = 0.1;
       //at 35 after (time>0): kp53F2 = 0.1;
       //rela KD T35
       //at 35 after (time>0): NFkB = 0.01;
       //at 35 after (time>0): kp65A =0.1;
       //at 35 after (time>0): kpNFkBB = 3;


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

    pycotools_mod = model.loada(TheModel.model_string, copasi_file=COPASI_FILE)
    pycotools_mod = tasks.TimeCourse(pycotools_mod, start=0, end=100).model


    if OPEN_WITH_COPASI:
        pycotools_mod.open()

    if PLOT_SIMULATION:
    
        NonSen = {
            'DDR': 0,
            'RAS': 0.5,
            'kDDRF': 0.01,
            'kRASF': 0.01,
            #'kDDRFRAS': 1,
        }
    
        plot_selection = {
            'Model C - Ras': ['RAS'],
            'Model C - Cell cycle_1': ['p21', 'p16','DDR'],
            'Model C - Cell cycle_2': ['p53', 'pp53', 'ppRb', 'pRb'],
            'Model C - Inflammation': ['p65IkBa', 'NFkB', 'pNFkB', 'IL6'],
            'Model C - Inflammation and p38': ['IL1a', 'CEBPb', 'p38', 'pp38'],
            'Model C - Notch signalling' : ['NOTCH1', 'NICD'],
        }

        ts = TimeSeries(mod.model_string, NonSen, 0, 100, 1000)
        tsp = TimeSeriesPlotter(ts, plot_selection, plot_dir=SIMULATION_DIRECTORY, fname='OIS_c_check.png',
                          savefig=True, legend=True, ncols=3, ylim=(-0.5,10))
        tsp.plot()
        tsp.data.to_csv('/Users/rebekahscanlan/python/tellurium/TGFmodel/data_files/OIS_c_check.csv')
        print(ts)
        
      
        
#        print(model_string.getSBML())
        
        
#        model = simplesbml.loadSBMLStr (mod.getSBML())
#        writeSBMLToFile(SBMLDocument, mod) 








