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
import random
from icecream import ic
import itertools
import json

def module(specs, new_tests: list):
    ic("module")
    ic(specs)
    ic(new_tests)
    new_list = []
    for new_test in new_tests:
        ic(new_test)
        for spec in specs:
            ic(spec)
            temp = new_test.copy()
            ic(temp)
            temp.update(spec)
            ic(temp)
            new_list.append(temp)
        ic(new_test)
        ic(new_tests)
    ic(new_list)
    
    return new_list
            
        

def test_thingo():
    a_list = [1, 2, 3, 4]
    b_list = a_list.copy()
    b_list.append([5, 6])
    ic(a_list)
    ic(b_list)
    
    first_test = dict()
    first_test = {"Category": "cat"}
    new_tests = [first_test]
    
    file_infos = [{"file_info": 1, "file_name": "something"}]
    new_tests = module(file_infos, new_tests)
    ic(new_tests)
    
    v_specs = [{"v_spec0": 123}, {"v_spec0": 234}]
    new_tests = module(v_specs, new_tests)
    ic(new_tests)
    
    ac = [{"ac": 10, "bg": 23}, {"ac": 20, "bg": 45}, {"ac": 30, "bg": 89}, {"ac": 40, "bg": 98}]
    new_tests = module(ac, new_tests)
    ic(new_tests)
    
    
    
if __name__ == '__main__':
    
    excel_input = "C:\grid_workspace\junction_rivers\src\junction_rivers\PSCAD\scripts\Config_test.xlsx"
    csv_output = "G:\Junction_Rivers\JRWF_PSCAD_Models\JRWF_PSCAD_SMIB_Siemens_GW_v5\spec_output.csv"
    mapping = "G:\Junction_Rivers\JRWF_PSCAD_Models\JRWF_PSCAD_SMIB_Siemens_GW_v5\jrwf_testbench_mapping.json"
    spec = SpecGenerator(excel_input, csv_output)
    
