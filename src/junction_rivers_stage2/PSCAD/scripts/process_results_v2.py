import logging
import os
import multiprocessing
import time
import argparse
import json
import pandas as pd
from rengen.plotting.Plotter import Plotter
from rengen.pscad import Psout
from rengen.spec.spec import load_specs_from_csv
from rengen.utils.time_utils import get_date_time_str
from icecream import ic

# from junction_rivers.Analysis import Analysis
from rengen.utils.gui_utils import prompt_for_multiple_filepaths, prompt_for_directory_path, std_script_title
from junction_rivers_stage2.plotting.JRWFStage2PlotterV1 import JRWFStage2Plotter

RUN_POST_PROCESSING = False

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("process_results_v2")
logger.setLevel(logging.INFO)


def per_scenario_postprocessing(spec_dict: dict, df: pd.DataFrame):
       
    # Per Scenario Analysis Added to JSON File
    spec_dict["analysis"] = {}
    
    return spec_dict, df
    
def process_results(
    src_data_path: os.PathLike,
    src_json_path: os.PathLike,
    pkl_path: os.PathLike,
    png_path: os.PathLike,
    pdf_path: os.PathLike,
    json_path: os.PathLike,
    secs_to_remove_from_traces: float,
    plotter: Plotter,
    delete_src_data: bool,
    run_post_processing: bool = RUN_POST_PROCESSING
):
    # Load Results Dataframe
    #logger.info(f"Reading: {src_data_path}")
    if ".psout" in src_data_path:
        psout = Psout(src_data_path)
        df = psout.to_df()
        
    elif ".pkl" in src_data_path:
        df = pd.read_pickle(src_data_path)
        
    else:
        logger.error("Invalid Source Data Extension")
        return  
    
    
    # Save Results to pkl file
    logger.info(f"Saving Data to: {pkl_path}")
    df[secs_to_remove_from_traces::].to_pickle(pkl_path)
    
    # Load Spec Dict from JSON 
    spec_dict = {}
    with open(src_json_path, 'r') as json_file:
        spec_dict = json.load(json_file)
    
    # # Run Postprocessing / Analysis
    if run_post_processing:
        #logger.info("Run Postprocessing...")
        try:
            spec_dict, df = per_scenario_postprocessing(spec_dict, df)
        except:
            logger.warning("per scenario postprocessing failed")
    
    # Run Plotter Script
    # logger.info("Plotting Results...")
    try:
        plotter.plot_from_df_and_dict(
            df=df,
            spec_dict=spec_dict,
            pdf_path=pdf_path,
            png_path=png_path,
        )
    except Exception as e:
        print(f"---------PLOTTER FAILED WITH EXCEPTION: {e}")
    
    #logger.info("Saving JSON Metadata...")
    with open(json_path, 'w') as json_file:
        json.dump(spec_dict, json_file, indent=2)
    
    if delete_src_data:
        #logger.info("Deleting Source Data...")
        os.remove(src_data_path)
        os.remove(src_json_path)
        
    return


def process_results_multi_thread(
    spec_path: os.PathLike,
    temp_results_dir: os.PathLike,
    results_dir: os.PathLike,
    processes: int = 1,
    secs_to_remove_from_traces: float = 0,
    delete_src_data: bool = False,
    data_source_extension: str = ".psout",
):
    plotter = JRWFStage2Plotter()

    spec = None
    if spec_path:
        while not os.path.exists(spec_path):
            time.sleep(30)
            logger.info(f"{spec_path} not found. Waiting 30 seconds...")

        logger.info(f"{spec_path} found.")
        ic(spec_path)
        spec = load_specs_from_csv(spec_path)
        ic(spec)

    # Results directory
    if spec_path is None:
       results_dir = os.path.join(results_dir, f"{get_date_time_str()}_results")
    print(f"\nResults Directory located at: \n   {results_dir}")
    os.makedirs(results_dir, exist_ok=True)

    # Intermetiate temp results directory
    # print(f"\nTemp Results Directory located at: \n   {temp_results_dir}")
    os.makedirs(temp_results_dir, exist_ok=True)

    # List of sudies if spec is provided
    study_list = []
    if spec_path:
        study_list = list(zip(list(spec["DIR"].values), list(spec["File_Name"].values)))
        num_studies = len(study_list)
    # While searching for results to process
    while True:
        #logger.info("Scanning for new results...")
        args = []
        # ic()
        # ic(args)
        for root, dirs, files in os.walk(temp_results_dir):
            relative_path = os.path.relpath(root, temp_results_dir)
            dst_results_dir = os.path.join(results_dir, relative_path)

            for file in files:
                file_ext_matched = os.path.splitext(file)[1] == data_source_extension
                file_base_name = os.path.splitext(file)[0]

                file_in_study_list = (relative_path, file_base_name) in study_list
                if file_ext_matched and (file_in_study_list or not spec_path):
                    src_data_path = os.path.join(root, file_base_name + data_source_extension)
                    src_json_path = os.path.join(root, file_base_name + ".json")

                    os.makedirs(dst_results_dir, exist_ok=True)

                    pkl_path = os.path.join(dst_results_dir, file_base_name + ".pkl")
                    png_path = os.path.join(dst_results_dir, file_base_name + ".png")
                    pdf_path = os.path.join(dst_results_dir, file_base_name + ".pdf")
                    json_path = os.path.join(dst_results_dir, file_base_name + ".json")

                    # Arguments for processing study
                    args.append([
                        src_data_path,
                        src_json_path,
                        pkl_path,
                        png_path,
                        pdf_path,
                        json_path,
                        secs_to_remove_from_traces,
                        plotter,
                        delete_src_data,
                    ])
                    # ic("________________")
                    for arg in args:
                        ic(arg[3])
                        # print("disabled ic")
                    # Remove study from study list if it was in there.
                    while (relative_path, file_base_name) in study_list:
                        study_list.remove((relative_path, file_base_name))
                    
        if args:
            remaining_str = f"({len(study_list)} of {num_studies} remaining)" if spec_path else ""
            logger.info(f"Processing {len(args)} Studies... {remaining_str}")

            # Multiprocess the processing of results
            if processes == 1:
                for proc_result_args in args:
                    process_results(*proc_result_args)
            else:
                with multiprocessing.Pool(processes=processes) as pool:
                    pool.starmap(process_results, args)

        # Stop searching after all studies removed from study list
        if study_list:
            # ic(study_list)
            time.sleep(30)
        else:
            logger.info("Processed all results")
            break
    plotter.plot_summary_pdf(results_dir)    
    


def process_results_single_thread(
    spec_path: os.PathLike,
    temp_results_dir: os.PathLike,
    results_dir: os.PathLike,
    secs_to_remove_from_traces: float = 0,
    delete_src_data: bool = False,
    data_source_extension: str = ".psout",
):
    process_results_multi_thread(spec_path, temp_results_dir, results_dir, 1, secs_to_remove_from_traces,
                                 delete_src_data, data_source_extension)
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', "--output-spec-path", type=str, default=None)
    parser.add_argument('-t', "--temp-results-dir", type=str, default=None)
    parser.add_argument('-r', "--results-dir", type=str, default=True)
    parser.add_argument('-p', "--processes", type=int, default=1, help="(number of processes for multiprocessing)")
    parser.add_argument('-k', "--secs-to-remove", type=float, default=0, help="(seconds removed from begining of stored pkl)")
    parser.add_argument('-d', "--delete-src-data", type=bool, default=False, help="(1 = delete source data, 0 = keep source data)")
    parser.add_argument('-e', "--data-source", type=str, default=".psout",help="(.out or .pkl)")
    args = parser.parse_args()

    logger.info("Started: process_results_v2.py")
    
    # TODO: the optional args have a minor bug, but to get it to work for now, use -t and -r and add the paths
    if args.results_dir is None:
        print("Please Select New Results Directory:")
        results_dir = prompt_for_directory_path(prompt_title="Select Results Dir.", initial_dir="G:\\Junction_Rivers\\JRWF_PSCAD_Models") #os.getcwd())
    else:
        results_dir = args.results_dir
    
    if args.temp_results_dir is None:
        print("Please Select Results Directory with PSOUT file:")
        temp_results_dir = prompt_for_directory_path(prompt_title="Select Results Dir. with PSOUT", initial_dir="D:\\Temp") #os.getcwd())
    else:
        temp_results_dir = args.temp_results_dir

    spec_path = args.output_spec_path    
    temp_results_dir = args.temp_results_dir
    results_dir = args.results_dir
    processes = args.processes
    secs_to_remove_from_traces = args.secs_to_remove
    delete_src_data = args.delete_src_data
    data_source_extension = args.data_source
    
    process_results_multi_thread(spec_path, temp_results_dir, results_dir, processes, secs_to_remove_from_traces,
                                 delete_src_data, data_source_extension)
    