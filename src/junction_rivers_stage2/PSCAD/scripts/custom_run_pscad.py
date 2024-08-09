from icecream import ic
import uuid
import pandas as pd
from rengen.pscad import launch_pscad, validate_pscad_model, PscadValidatorBehaviour, run_pscad_spec
from rengen.spec.spec import load_specs_from_csv
# from rengen.plotting.UbwfAPscadPlotterV1 import UbwfAPscadPlotterV1
from rengen.utils.gui_utils import prompt_for_multiple_filepaths, prompt_for_directory_path
from junction_rivers_stage2.PSCAD.scripts.process_results_v2 import process_results_single_thread
from junction_rivers_stage2.plotting.JRWFStage2PlotterV1 import JRWFStage2Plotter
import os
from datetime import datetime

from multiprocessing import Pool
import tempfile
import json

def make_temp_directory(parent_dir) -> os.PathLike:
    unique_id = str(uuid.uuid4())[:32]
    tmp_path = os.path.join(parent_dir, unique_id)
    os.makedirs(tmp_path)
    return tmp_path


def get_date_time_str() -> str:
    now = datetime.now()
    # return now.strftime("%Y_%m_%d_%H%M")
    return now.strftime("%Y_%m_%d")

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


def launch_pscad_study(
    study_title: str,
    model_location: os.PathLike,
    spec_files: os.PathLike = None,
    results_dir: os.PathLike = r"G:\Junction_Rivers_Stage2\JRWF_S2_PSCAD_Models\JRWF_S2_PSCAD_GW_EKS_v0.10-SP\results_gs",
    model_tmp_dir: os.PathLike = "D:/Temp/",
    volley_size: int = 13,
    group: int = None,
    file_index: int = None,
    gs_overrides: dict = None
):
    
    """-------------------------------------------------------------------------------  
       Function To Run PSCAD Simulation Using the Regen Codebase
    -------------------------------------------------------------------------------"""

    print(
        study_title,
        model_location,
        spec_files,
        "-----------------------------------------------",
    )



    print(model_location)

    # Derived Parameters
    workspace_path = os.path.join(model_location,"model_gs_test/JRWF_S2_SMIB_Workspace.pswx")
    tuning_file = os.path.join(model_location,"Global_Substitutions/jrwf_s2_gs.csv")
    mapping_file = os.path.join(model_location,"jrwf_testbench_mapping.json")

    print(workspace_path)

    # Propt for Spec Files
    if spec_files is None:
        spec_files = prompt_for_multiple_filepaths(prompt_title="Select Spec Files.")
        print("The following specs have been selected to run:")
        for file in spec_files:
            print(file)

    # Load Spec Files
    spec = load_specs_from_csv(spec_files)

    # Propt for Results Directory
    if results_dir is None:
        # print("Please Select Results Directory.")
        # results_dir = prompt_for_directory_path("Select Output/Results Directory")
        results_dir = os.path.join(model_location,f"../results/{get_date_time_str()}_{study_title}_Results")

    print(f"Results Directory located at: {results_dir}")

    # Make Temp Directory
    temp_results_dir = make_temp_directory(model_tmp_dir)
    print(f"Model Temp Directory located at: {temp_results_dir}")

    # Filter Specs by Group and Index
    filtered_spec = spec
    if group is not None:
        filtered_spec = spec[spec['Grouping'] == group]
    if file_index is not None:
        str_list = file_index.split(",")
        suffixes = [x.zfill(3) for x in str_list]
        filtered_spec = filtered_spec[filtered_spec['File_Name'].str.endswith(tuple(suffixes))]
    
    output_spec_title = "output_spec.csv"
    output_spec_path = os.path.join(results_dir, output_spec_title)
    print(f"Setting up {len(filtered_spec)} simulations:")

    # Lauch PSCAD Instance
    pscad, project = launch_pscad(workspace_path=workspace_path, project_name="JRWF_S2_SMIB", copy_to_dir=temp_results_dir)

    tuning_df = dataframe_from_gs_csv(tuning_file)
    tuning_dict = tuning_df.set_index('Name ')[' Value '].to_dict()
    
    
    
    with tempfile.NamedTemporaryFile(delete=False, mode='w+') as temp_file:
        # Get the name of the temporary file
        temp_json_tuning_file = temp_file.name
        print(f"Temporary file created: {temp_json_tuning_file}")
        json.dump(tuning_dict, temp_file, indent=4)
        temp_file.flush()

    # Overwrite PSCAD Global Subs with Tuning File
    # validate_pscad_model(project=project, json_path=tuning_file, behaviour=PscadValidatorBehaviour.OVERWRITE_PSCAD)
    validate_pscad_model(project=project, json_path=temp_json_tuning_file, behaviour=PscadValidatorBehaviour.OVERWRITE_PSCAD)
    
    # Launch Study
    run_pscad_spec(
        pscad=pscad,
        project=project,
        spec=filtered_spec,
        testbench_map_path=mapping_file,
        gs_overrides=gs_overrides,
        plotter = JRWFStage2Plotter,
        volley_size=volley_size, 
        results_dir=temp_results_dir,
        allow_existing_results_dir=True,
        output_spec_path=output_spec_path,
        prioritised_column=' Category ',
        quit_after=False,
    )
    # run_pscad_spec(
    #     pscad=pscad,
    #     project=project,
    #     spec=filtered_spec,
    #     testbench_map_path=mapping_file,
    #     volley_size=volley_size,  # Empirically determined, Project/Model/PC Specifric.
    #     plotter=UbwfAPscadPlotterV1(),
    #     quit_after=True,
    #     results_dir=results_dir,
    #     allow_existing_results_dir=True,
    #     subfolder_column='DIR',
    #     filename_column='File_Name',
    #     output_spec_path="output_spec.csv",
    #     save_output_as_pickle=True,
    #     prioritised_column='Category', 
    #     prioritised_values='Fault_Ohms',
    #     gs_overrides=gs_overrides,
    # )

def function_wrapper(par_dict):
    launch_pscad_study(**par_dict)

if __name__ == '__main__':


    parameter_list = []

    parameter_list.append(
        {
            "study_title": f"Case_1",
            "model_location": "G:\Junction_Rivers_Stage2\JRWF_S2_PSCAD_Models\JRWF_S2_PSCAD_GW_EKS_v0.10-SP",
            "spec_files": r"G:\Junction_Rivers_Stage2\JRWF_S2_PSCAD_Models\JRWF_S2_PSCAD_GW_EKS_v0.10-SP\specs\2024_08_06_1433_runner_spec_flatruns_gs.csv",
            # "group": 1,
            r"gs_overrides": {
                "WT_Ki_D21": 10,
                "WT_QVMode_D22": 2,
                "WT_Kp_D23": 0.5,
                "model_v": 0.1001,
            }
        },
    )


    parameter_list.append(
        {
            "study_title": f"Case_1",
            "model_location": "G:\Junction_Rivers_Stage2\JRWF_S2_PSCAD_Models\JRWF_S2_PSCAD_GW_EKS_v0.10-SP",
            "spec_files": r"G:\Junction_Rivers_Stage2\JRWF_S2_PSCAD_Models\JRWF_S2_PSCAD_GW_EKS_v0.10-SP\specs\2024_08_06_1433_runner_spec_flatruns_gs.csv",
            # "group": 1,
            r"gs_overrides": {
                "WT_Ki_D21": 10,
                "WT_QVMode_D22": 2,
                "WT_Kp_D23": 1,
                "model_v": 0.1001,
            }
        },
    )

    # parameter_list.append(
    #     {
    #         "study_title": f"Case_1",
    #         "model_location": "G:\Junction_Rivers_Stage2\JRWF_S2_PSCAD_Models\JRWF_S2_PSCAD_GW_EKS_v0.10-SP",
    #         "spec_files": r"G:\Junction_Rivers_Stage2\JRWF_S2_PSCAD_Models\JRWF_S2_PSCAD_GW_EKS_v0.10-SP\specs\2024_08_06_1433_runner_spec_flatruns_gs.csv",
    #         # "group": 1,
    #         r"gs_overrides": {
    #             "WT_Ki_D21": 10,
    #             "WT_QVMode_D22": 2,
    #             "WT_Kp_D23": 1.5,
    #             "model_v": 0.1001,
    #         }
    #     },
    # )


    # parameter_list.append(
    #     {
    #         "study_title": f"Case_1",
    #         "model_location": "G:\Junction_Rivers_Stage2\JRWF_S2_PSCAD_Models\JRWF_S2_PSCAD_GW_EKS_v0.10-SP",
    #         "spec_files": r"G:\Junction_Rivers_Stage2\JRWF_S2_PSCAD_Models\JRWF_S2_PSCAD_GW_EKS_v0.10-SP\specs\2024_08_06_1433_runner_spec_flatruns.csv",
    #         # "group": 1,
    #         r"gs_overrides": {
    #             "WT_Ki_D21": 10,
    #             "WT_QVMode_D22": 2,
    #             "WT_Kp_D23": 2,
    #         }
    #     },
    # )


    # parameter_list.append(
    #     {
    #         "study_title": f"Case_2", 
    #         "model_location": "D:/grid_workspace/tasks/Operational_Limits_Study_2/Models/GBWF-v3.25-Case_2",
    #         "spec_files": r"D:/grid_workspace/tasks/Operational_Limits_Study_2/specs/spec__DMAT_3215_OLS_Case_2.csv",
    #         # "group": 2,
    #         "gs_overrides": {
    #             "BESS1_PCS_Inverters_Online": 0,
    #             "BESS2_PCS_Inverters_Online": 16,
    #             "BESS3_PCS_Inverters_Online": 0,
    #             "BESS4_PCS_Inverters_Online": 0,
    #             "WTG1_Units_Online": 0,
    #             "WTG2_Units_Online": 19,
    #             "WTG3_Units_Online": 0,
    #             "WTG4_Units_Online": 0,
    #             }
    #     },
    # )

    # parameter_list.append(
    #     {
    #         "study_title": "Case_3", 
    #         "model_location": "D:/grid_workspace/tasks/Operational_Limits_Study_2/Models/GBWF-v3.25-Case_3",
    #         "spec_files": r"D:/grid_workspace/tasks/Operational_Limits_Study_2/specs/spec__DMAT_3215_OLS_Case_3.csv",
    #         # "group": 2,
    #         "gs_overrides": {
    #             "BESS1_PCS_Inverters_Online": 32,
    #             "BESS2_PCS_Inverters_Online": 16,
    #             "BESS3_PCS_Inverters_Online": 32,
    #             "BESS4_PCS_Inverters_Online": 32
    #         }
    #     },
    # )


    with Pool(3) as pool:
        pool.map(function_wrapper, parameter_list)
