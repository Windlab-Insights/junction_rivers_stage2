import multiprocessing
from functools import partial
import shutil
import subprocess
import time
from multiprocessing import Pool
from junction_rivers.PSCAD.scripts.process_results import process_results_single_thread
from icecream import ic
import argparse
import json
import pandas as pd
from rengen.pscad import Psout, launch_pscad, run_pscad_spec, validate_pscad_model, PscadValidatorBehaviour
from rengen.spec.spec import load_specs_from_csv
from rengen.utils.gui_utils import prompt_for_multiple_filepaths, prompt_for_directory_path, std_script_title
import os
import tempfile

from rengen.utils.os_utils import make_temp_directory, file_exists_beneath_directory
from rengen.utils.time_utils import get_date_time_str

MODEL_PATH = "."
DEFAULT_VOLLEY = 12
PROJECT_NAME = "JRWF_SMIB"


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


def dataframe_from_gs_csv(tuning_file):
    # Function to determine if a value is a number
    def convert_type(value):
        if isinstance(value, str):

            # Try to convert to integer
            try:
                return int(value)
            except ValueError:
                pass

            # Try to convert to float
            try:
                return float(value)
            except ValueError:
                pass

        # if already a string, or converting to int/float failed, return the original value.
        return value

    # Read the CSV file into a DataFrame with the first row as header
    df = pd.read_csv(tuning_file, header=0)

    # Apply the type conversion function to the ' Value ' column
    df[' Value '] = df[' Value '].apply(convert_type)
    return df


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', "--group", type=int, default=None)
    parser.add_argument('-v', "--volley-size", type=int, default=DEFAULT_VOLLEY)
    parser.add_argument('-s', "--spec-path", type=str, default=None)
    parser.add_argument('-t', "--temp-dir", type=str, default="C:\\Temp")
    parser.add_argument('-r', "--results-dir", type=str, default=None)
    parser.add_argument('-f', "--file-index", type=str, default=None, help="Accepts comma-separated list. Numbers correspond to excel index.")
    parser.add_argument('-w', '--work-partition', action=WorkerValidation, dest='worker_partition', help='Partitions Work. Use pattern \"worker_total:worker_id\". Where worker_total=2,3,4,5, and worker_id=1,..,worker_total')
    parser.add_argument('-m', "--missing-only", action='store_true', default=False, help="Given a results directory, only runs tests that do not currently exist.")
    parser.add_argument("--store-init-traces", action='store_true', default=True, help="If set, the init. period will NOT be removed from the output traces.")
    parser.add_argument("--minimum-traces", action='store_true', default=False, help="If set, the only the minimum traces needed for the grid-plot will be stored in the pkl file.")
    parser.add_argument("-p", '--model-path', default="G:\Junction_Rivers\JRWF_PSCAD_Models\JRWF_PSCAD_SMIB_Siemens_GW_v5")
    args = parser.parse_args()

    ic(args)

    volley_size = args.volley_size
    
    if args.model_path is None:
        print("Please Select Desired OEM Model Path To Run:")
        model_path = prompt_for_multiple_filepaths(prompt_title="Select Model Path", initial_dir="G:\Junction_Rivers\JRWF_PSCAD_Models")
    else: 
        model_path = args.model_path
    
    workspace_path = os.path.join(model_path, "model/JRWF_SMIB_Workspace.pswx")
    tuning_file = os.path.join(model_path,"Global_Substitutions/jrwf_gs.csv")
    mapping_file = os.path.join(model_path,"jrwf_testbench_mapping.json")

    std_script_title("run_pscad_simulation_script.py")

    file_paths = []
    if args.spec_path is None:
        print("Please Select Specs To Run:")
        file_paths = prompt_for_multiple_filepaths(prompt_title="Select Spec Files.", initial_dir=model_path) #os.getcwd())
    else:
        file_paths = [args.spec_path]
    for file_path in file_paths:
            print(f"   {file_path}")
    spec = load_specs_from_csv(file_paths)

    if args.temp_dir is None:
        print("\nPlease Select Temp Directory (Not within dropbox):")
        model_temp_dir = prompt_for_directory_path("Select Temp. Directory", initial_dir=os.getcwd())
    else:
        model_temp_dir = args.temp_dir
    temp_path = make_temp_directory(model_temp_dir)
    print(f"\nModel Temp Directory located at: \n   {temp_path}")

    if args.results_dir is None:
        print("\nPlease Select Results Directory.")
        results_base_dir = prompt_for_directory_path("Select Output/Results Directory", initial_dir=model_path)
    else:
        results_base_dir = args.results_dir
    # os.makedirs(results_base_dir)
    
    results_dir = os.path.join(results_base_dir,f"{get_date_time_str()}_results")
    print(f"\nResults Directory located at: \n   {results_dir}")
    os.makedirs(results_dir,exist_ok=True)

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

    temp_results_dir = os.path.join(temp_path, "__temp_results")
    print(f"\nTemp Results Directory located at: \n   {temp_results_dir}")
    
    output_spec_path = os.path.join(results_dir, "output_spec.csv")

    pool = multiprocessing.Pool(processes=2)
    pool.apply_async(process_results_single_thread, [output_spec_path, temp_results_dir, results_dir])

    print(f"Workspace_Path:{workspace_path}")
    pscad, project = launch_pscad(workspace_path=workspace_path, project_name=PROJECT_NAME, copy_to_dir=temp_path)
    tuning_df = dataframe_from_gs_csv(tuning_file)
    tuning_dict = tuning_df.set_index('Name ')[' Value '].to_dict()
    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False, mode='w+') as temp_file:
        # Get the name of the temporary file
        temp_json_tuning_file = temp_file.name
        print(f"Temporary file created: {temp_json_tuning_file}")
        json.dump(tuning_dict, temp_file, indent=4)
        temp_file.flush()

    validate_pscad_model(project=project, json_path=temp_json_tuning_file, behaviour=PscadValidatorBehaviour.OVERWRITE_PSCAD)
    os.remove(temp_json_tuning_file)

    run_pscad_spec(
        pscad=pscad,
        project=project,
        spec=filtered_spec,
        testbench_map_path=mapping_file,
        volley_size=volley_size, 
        results_dir=temp_results_dir,
        allow_existing_results_dir=True,
        output_spec_path=output_spec_path,
        prioritised_column=' Category ',
        quit_after=False,
    )
    
    print("Finished Run.")
    
    pool.close()
    pool.join()
    