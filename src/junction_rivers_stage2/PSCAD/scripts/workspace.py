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
    
    excel_input = "C:\grid_workspace\junction_rivers_stage2\src\junction_rivers_stage2\PSCAD\scripts\Config_test.xlsx"
    csv_output_EKS = "G:\Junction_Rivers_Stage2\JRWF_S2_PSCAD_Models\JRWF_S2_PSCAD_GW_EKS_v0.9\specs\\runner_spec.csv"
    spec = SpecGenerator(excel_input, csv_output_EKS)
    
