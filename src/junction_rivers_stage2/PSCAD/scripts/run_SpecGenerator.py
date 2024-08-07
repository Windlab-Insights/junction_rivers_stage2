from collections import OrderedDict
import pandas as pd
from junction_rivers.Analysis import Analysis
import os
import argparse
import numpy as np
from typing import List, Optional, Tuple, Union, Dict
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as backend_pdf
import matplotlib.gridspec as gridspec
from SpecGenerator import SpecGenerator
from rengen.utils.gui_utils import prompt_for_directory_path,prompt_for_open_filepath
import math
import random
from icecream import ic
import itertools
import json   
from rengen.utils.time_utils import get_date_time_str
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-suite', '--suite', choices=['dmat', 'r_dmat', 'csr', 'none'], type=str)
    parser.add_argument('-all', action='store_true', help='Run all tests regardless of action')
    parser.add_argument('-cust', action='store_true', help='Run only the tests with action')
    parser.add_argument('-c', "--config_test_path", type=str, default="C:\grid_workspace\junction_rivers_stage2\src\junction_rivers_stage2\PSCAD\scripts\Config_test.xlsx")
    parser.add_argument('-s', "--spec_output_direc", type=str, default="G:\Junction_Rivers_Stage2\JRWF_S2_PSCAD_Models\JRWF_S2_PSCAD_GW_EKS_v0.3\specs")
    
    args = parser.parse_args()
        
    config_test = []
    if args.config_test_path is None:
        print("Please Select the Config_test.xlsx file to run:")
        config_test = prompt_for_open_filepath(prompt_title="Select Config_test.xlsx") #os.getcwd())
    else:
        config_test = args.config_test_path
    
    # spec_output_dir = []
    initial_spec_path = "G:\Junction_Rivers_Stage2\JRWF_S2_PSCAD_Models"
    if args.spec_output_direc is None:
        print("Please Select the directory for runner_spec.csv to be saved:")
        spec_output_dir = prompt_for_directory_path(prompt_title="Select runner_spec.csv directory:",initial_dir=initial_spec_path) #os.getcwd())
    else:
        spec_output_dir = args.spec_output_direc
    spec_output_path = os.path.join(spec_output_dir,f"{get_date_time_str()}_runner_spec.csv")
    
    if args.suite is None:
        suite = None
        while not (suite == "dmat" or suite == "r_dmat" or suite == "csr" or suite == "none"):
            suite = input("Please specify suite (dmat, r_dmat, csr, none): ")
    else:
        suite = args.suite
    
    if args.all:
        full_suite = True
    elif args.cust:
        full_suite = False
    else:
        run_all = None
        while not (run_all == "yes" or run_all == "no"):
            run_all = input("Would you like to run all the tests regardless of action? (yes, no): ")
        full_suite = (run_all == "yes")
        
    ic(suite)
    ic(full_suite)
    spec = SpecGenerator(config_test, spec_output_dir, suite, full_suite)
    
