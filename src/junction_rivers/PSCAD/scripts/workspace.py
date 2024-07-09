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
    
if __name__ == '__main__':
    
    excel_input = "C:\grid_workspace\junction_rivers\src\junction_rivers\PSCAD\scripts\Config_test.xlsx"
    csv_output_Siemens = "G:\Junction_Rivers\JRWF_PSCAD_Models\JRWF_PSCAD_SMIB_Siemens_GW_v5\specs\\runner_spec.csv"
    csv_output_EKS = "G:\Junction_Rivers\JRWF_PSCAD_Models\JRWF_PSCAD_SMIB_EKS_GW\specs\\runner_spec.csv"
    mapping = "G:\Junction_Rivers\JRWF_PSCAD_Models\jrwf_testbench_mapping.json"
    spec = SpecGenerator(excel_input, csv_output_Siemens)
    
