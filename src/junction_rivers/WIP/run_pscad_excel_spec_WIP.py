import csv
import multiprocessing
from functools import partial
import shutil
import subprocess
import time
from multiprocessing import Pool
from icecream import ic
import argparse
import json
import pandas as pd
from rengen.pscad import Psout, launch_pscad, run_pscad_spec
from rengen.spec.spec import load_specs_from_csv
from rengen.plotting.BBWFPlotterV2 import BBWFPlotter
from rengen.utils.gui_utils import prompt_for_multiple_filepaths, prompt_for_directory_path, std_script_title
import os

import openpyxl as opx

from rengen.utils.helper_functions import get_dict_from_excel_table
from rengen.utils.os_utils import make_temp_directory, file_exists_beneath_directory
from rengen.utils.time_utils import get_date_time_str

MODEL_PATH = "."
DEFAULT_VOLLEY = 10

# Custom action class for validation
class WorkerValidation(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # Perform your validation logic here
        input = values.split(":")
        worker_total = int(input[0])
        worker_id = int(input[1])
        if worker_total not in [2, 3, 4, 5]:
            parser.error(f"The worker_total must be 2,3,4,5: Was provided with {worker_total}")
        if worker_id < 1 or worker_id > worker_total:
            parser.error(f"The workerid must be between 1, 2, .., worker_total: was provided with {worker_id}")
        setattr(namespace, self.dest, {'worker_total':worker_total, 'worker_id':worker_id})
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', "--group", type=int, default=None)
    parser.add_argument('-v', "--volley-size", type=int, default=DEFAULT_VOLLEY)
    parser.add_argument('-e', "--excel-spec-path", type=str, default=None, help="Excel spec file path (ie. \"./spec.xlsx\" )")
    parser.add_argument('-s', "--selected-studies", type=str, default=None, help="selected studies (ie \"csr_s5255\")")
    parser.add_argument('-t', "--temp-dir", type=str, default="D:\\Temp")
    parser.add_argument('-r', "--results-dir", type=str, default=".\\results")
    parser.add_argument('-f', "--file-index", type=str, default=None, help="Accepts comma-separated list. Numbers correspond to excel index.")
    parser.add_argument('-w', '--work-partition', action=WorkerValidation, dest='worker_partition', help='Partitions Work. Use pattern \"worker_total:worker_id\". Where worker_total=2,3,4,5, and worker_id=1,..,worker_total')
    parser.add_argument('-m', "--missing-only", action='store_true', default=False, help="Given a results directory, only runs tests that do not currently exist.")
    parser.add_argument("--store-init-traces", action='store_true', default=True, help="If set, the init. period will NOT be removed from the output traces.")
    parser.add_argument("--minimum-traces", action='store_true', default=False, help="If set, the only the minimum traces needed for the grid-plot will be stored in the pkl file.")
    args = parser.parse_args()

    ic(args)

    volley_size = args.volley_size
    excel_spec_path = args.excel_spec_path
    selected_studies = str(args.selected_studies).strip(" []").split(",")

    workspace_path = os.path.join(MODEL_PATH, "model/BBWF_WS.pswx")
    tuning_file = os.path.join(MODEL_PATH,"bbwf_tuning.json")
    mapping_file = os.path.join(MODEL_PATH,"bbwf_testbench_mapping.json")

    std_script_title("run_pscad_simulation_script.py")


    if args.temp_dir is None:
        print("\nPlease Select Temp Directory (Not within dropbox):")
        model_temp_dir = prompt_for_directory_path("Select Temp. Directory", initial_dir=os.getcwd())
    else:
        model_temp_dir = args.temp_dir
    temp_path = make_temp_directory(model_temp_dir)
    print(f"\nModel Temp Directory located at: \n   {temp_path}")


    if args.results_dir is None:
        print("\nPlease Select Results Directory.")
        results_base_dir = prompt_for_directory_path("Select Output/Results Directory", initial_dir=os.getcwd())
    else:
        results_base_dir = args.results_dir
        
    results_dir = os.path.join(results_base_dir,f"{get_date_time_str()}_results")
    print(f"\nResults Directory located at: \n   {results_dir}")
    os.makedirs(results_dir,exist_ok=True)
    
    
    
    spec_wb = opx.load_workbook(excel_spec_path)
    
    
    file_paths = []
    
    for selected_study in selected_studies:
        
        file_path = f"{temp_path}/{selected_study}.csv"
        with open(file_path, 'w') as file:
            csv_writer = csv.writer(file)
            
            ws =spec_wb[selected_study]
        
            for row in ws.rows:
                
                csv_writer.writerow([cell.value for cell in row])
            
            file_paths.append(file_path)

    spec = load_specs_from_csv(file_paths)
    

    filtered_spec = spec

    if args.worker_partition is not None:
        worker_total = args.worker_partition['worker_total']
        worker_id = args.worker_partition['worker_id']

        worker_id_column = f"_{worker_total}WorkerPartition"
        if worker_id_column not in filtered_spec.columns:
            print(f"Spec needs {worker_id_column} column, in order to divide work amongst workers.")
            exit()
        else:
            filtered_spec = filtered_spec[filtered_spec[worker_id_column] == worker_id]

    if args.group is not None:
        filtered_spec = spec[spec['Grouping'] == args.group]

    if args.file_index is not None:
        str_list = args.file_index.split(",")
        suffixes = [x.zfill(3) for x in str_list]
        filtered_spec = filtered_spec[filtered_spec['File_Name'].str.endswith(tuple(suffixes))]

    if args.missing_only:
        
        print("Please Select Existing Results Root Directory.")
        results_root_dir = prompt_for_directory_path("Select Root Results Directory", initial_dir=os.getcwd())
        
        filtered_spec['_File_Exists'] = False
        for index, row in filtered_spec.iterrows():
            pkl_exists = file_exists_beneath_directory(results_root_dir, f"{row['File_Name']}.pkl")
            psout_exists = file_exists_beneath_directory(results_root_dir, f"{row['File_Name']}.psout")
            if pkl_exists or psout_exists:
                filtered_spec.at[index, '_File_Exists'] = True
        filtered_spec = filtered_spec[filtered_spec['_File_Exists'] == False]

    print(f"\nSetting up {len(filtered_spec)} simulations.")

    temp_results_dir = os.path.join(temp_path,"__temp_results")
    print(f"\nTemp Results Directory located at: \n   {temp_results_dir}")
    
    output_spec_path = os.path.join(results_dir,"output_spec.csv")


    pool = multiprocessing.Pool(processes=2) 
    
    pool.apply_async(subprocess.run, [["py","./process_results.py",'-s',output_spec_path,'-t',temp_results_dir,'-r',results_dir]])
    
    pscad, project = launch_pscad(workspace_path=workspace_path, project_name="BBWF_Project", copy_to_dir=temp_path)
    
    run_pscad_spec(
        pscad=pscad,
        project=project,
        spec=filtered_spec,
        testbench_map_path=mapping_file,
        volley_size=volley_size, 
        results_dir=temp_results_dir,
        allow_existing_results_dir=True,
        subfolder_column='DIR',
        filename_column='File_Name',
        output_spec_path=output_spec_path,
        prioritised_column='Category',
    )
    
    pool.close()
    
    pool.join()
    