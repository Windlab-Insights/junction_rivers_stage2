from collections import OrderedDict
import pandas as pd
from junction_rivers.Analysis import Analysis
import os
import numpy as np
from typing import List, Optional, Tuple, Union, Dict
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as backend_pdf
import matplotlib.gridspec as gridspec
from SpecGenerator import SpecGenerator
import math

def calc_vslack_from_vpoc():
    vpoc = 0.9
    s_nom = 100
    ppoc = 150/s_nom
    qpoc = 59/s_nom
    spoc = math.sqrt(ppoc**2 + qpoc**2)
    zbase = 220**2/s_nom
    
    v_nom = 220
    (rgrid, xgrid, zgrid) = (10.3758, 0.1453 * 2*50 * math.pi, math.sqrt(10.3758**2 + (0.1453 * 2*50 * math.pi)**2))
    rgrid_pu = rgrid/ zbase
    xgrid_pu = xgrid/ zbase
    zgrid_pu = zgrid/ zbase
    print(f'### rgrid_pu = {rgrid_pu}, xgrid_pu = {xgrid_pu}, zgrid_pu = {zgrid_pu}')
    vslack = math.sqrt(vpoc**2 + spoc**2/vpoc**2*zgrid_pu**2 - 2*(ppoc*rgrid_pu + qpoc*xgrid_pu))
    return vslack

def calc_vpoc_from_vslack():
    vslack = 0.81511
    s_nom = 100
    ppoc = 150/s_nom
    qpoc = 59/s_nom
    spoc = math.sqrt(ppoc**2 + qpoc**2)
    zbase = 220**2/s_nom
    (rgrid, xgrid, zgrid) = (10.3758, 0.1453 * 2*50 * math.pi, math.sqrt(10.3758**2 + (0.1453 * 2*50 * math.pi)**2))
    rgrid_pu = rgrid/ zbase
    xgrid_pu = xgrid/ zbase
    zgrid_pu = zgrid/ zbase
    print(f'### rgrid_pu = {rgrid_pu}, xgrid_pu = {xgrid_pu}, zgrid_pu = {zgrid_pu}')
    vpoc = math.sqrt((vslack**2 + 2*(ppoc*rgrid_pu + qpoc*xgrid_pu) + math.sqrt((vslack**2 + 2*(ppoc*rgrid_pu + qpoc*xgrid_pu))**2 - 4*spoc**2*zgrid_pu**2))/2)
    print(f'### vgrid = {vpoc}')
    return vpoc
    
if __name__ == '__main__':
    # vslack = calc_vslack_from_vpoc()
    # # vslack =0
    # vpoc = calc_vpoc_from_vslack()
    # # vpoc = 0
    # print(f"vslack = {vslack}, vpoc = {vpoc}")
    
    excel_input = "C:\grid_workspace\junction_rivers\src\junction_rivers\PSCAD\scripts\Config_test.xlsx"
    csv_output = "C:\Temp\csv_output.csv"
    mapping = "G:\Junction_Rivers\JRWF_PSCAD_Models\JRWF_PSCAD_SMIB_Siemens_GW_v5\jrwf_testbench_mapping.json"
    spec = SpecGenerator(excel_input, csv_output)
    
    
##################### CALCULATIONS #######################
        
    # def calc_fault_level(self, row):
    #     scr = row['SCR']
    #     print(f'### scr = {scr}')
    #     fault_level = scr * self.p_nom
    #     print(f'### fault_level = {fault_level}')
    #     return fault_level
    
    # def calc_grid_impedence(self, row):
    #     scr = row['SCR']
    #     x2r = row['X2R']
    #     print(f'### scr = {scr}, x2r = {x2r}')
    #     fl = self.calc_fault_level(scr)
    #     r = self.v_nom**2 * fl / math.sqrt(1+x2r**2)
    #     x = r * x2r
    #     z = math.sqrt(r**2 + x**2)
    #     print(f'### r = {r}, x = {x}, z = {z}')
    #     return (r, x, z)
        