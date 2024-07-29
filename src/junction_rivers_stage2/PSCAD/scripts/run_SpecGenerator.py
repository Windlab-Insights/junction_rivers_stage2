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
    parser.add_argument("-p", '--model-path', type=str, default=None)
    parser.add_argument('-c', "--config_test_path", type=str, default=None)
    parser.add_argument('-s', "--spec_output_direc", type=str, default=None)
    
    args = parser.parse_args()
    
    if args.model_path is None:
        print("Please Select Model Path:")
        model_path = prompt_for_directory_path(prompt_title="Select Model Path", initial_dir="G:\Junction_Rivers\JRWF_PSCAD_Models")
    else: 
        model_path = args.model_path
    
    config_test = []
    if args.config_test_path is None:
        print("Please Select the Config_test.xlsx file to run:")
        config_test = prompt_for_open_filepath(prompt_title="Select Config_test.xlsx") #os.getcwd())
    else:
        config_test = args.config_test_path
    
    # spec_output_dir = []
    if args.spec_output_direc is None:
        print("Please Select the directory for runner_spec.csv to be saved:")
        spec_output_dir = prompt_for_directory_path(prompt_title="Select runner_spec.csv directory:",initial_dir=model_path) #os.getcwd())
    else:
        spec_output_dir = args.spec_output_direc
    spec_output_path = os.path.join(spec_output_dir,f"{get_date_time_str()}_runner_spec.csv")
    
    spec = SpecGenerator(config_test, spec_output_path)
    
