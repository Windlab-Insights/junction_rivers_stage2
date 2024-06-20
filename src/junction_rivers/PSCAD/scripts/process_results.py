import logging
import os
import multiprocessing
import time
import argparse
import json
import pandas as pd
from collections import OrderedDict
from rengen.plotting.Plotter import Plotter
from rengen.pscad import Psout
from rengen.spec.spec import load_specs_from_csv
from rengen.utils.time_utils import get_date_time_str
from PyPDF2 import PdfWriter, PdfReader

from junction_rivers.analysis.analysis__csr_s5255_low_voltage_faults import per_scenario_analysis as s5255_low_voltage_fault_analysis
from junction_rivers.analysis.analysis__csr_s5255_high_voltage_faults import per_scenario_analysis as s5255_high_voltage_fault_analysis
from junction_rivers.analysis.analysis__csr_s52514_pref_step import s52514_pref_step_analysis
from rengen.utils.gui_utils import prompt_for_multiple_filepaths, prompt_for_directory_path, std_script_title
from junction_rivers.analysis.analysis__csr_s5253_s5258_freq_dist import s5253_s5258_freq_dist_analysis
from junction_rivers.analysis.signal_analysis import get_expected_fdroop_signal, get_expected_vdroop_signal
from junction_rivers.Analysis import Analysis

from junction_rivers.plotting.JRWFPlotterV2 import JRWFPlotter

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("process_results")
logger.setLevel(logging.INFO)


def per_scenario_postprocessing(spec_dict: dict, df: pd.DataFrame):
    
    # Signals Added to Dataframe
    # df["POC_Pref_Fdroop_MW"] = get_expected_fdroop_signal(spec_dict, df)
    # df["Qref_droop_Ideal_MVAr"]= get_expected_vdroop_signal(spec_dict, df)
       
    # Per Scenario Analysis Added to JSON File
    spec_dict["analysis"] = {}
        
    if "s5255_Fault_Ures" in spec_dict["File_Name"] or "s5255_Fault_Ohms" in spec_dict["File_Name"]:
        
            spec_dict["analysis"].update(
                s5255_low_voltage_fault_analysis(spec_dict["substitutions"], spec_dict, df)
            )
            
    if "s5255_Fault_TOV" in spec_dict["File_Name"]:
        
        spec_dict["analysis"].update(
            s5255_high_voltage_fault_analysis(spec_dict["substitutions"], spec_dict, df)
        )
        
    if "s5258" in spec_dict["File_Name"] or "s5253" in spec_dict["File_Name"]:
        
        spec_dict["analysis"].update(
            s5253_s5258_freq_dist_analysis(spec_dict, df)
        )
    
    if "Pref" in spec_dict["File_Name"]:
        
        spec_dict["analysis"].update(
            s52514_pref_step_analysis(spec_dict, df)
        )
    
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
    analysis: Analysis,
    delete_src_data: bool,
):
    
    # Load Results Dataframe
    logger.info(f"Reading: {src_data_path}")
    if ".psout" in src_data_path:
        try:
            psout = Psout(src_data_path)
            df = psout.to_df()
        except Exception as e:
            print(f"### psout read failed {e}")
        
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
    
    # Run Postprocessing / Analysis
    try:
        spec_dict, df = per_scenario_postprocessing(spec_dict, df)
    except:
        logger.warning("per scenario postprocessing failed")
    
    # Run Plotter Script
    logger.info("Plotting Results...")
    print("### plotting started consider:" + pdf_path)
    try:
        plotter.plot_from_df_and_dict(
            df=df,
            spec_dict=spec_dict,
            pdf_path=pdf_path,
            png_path=png_path,
        )
        print("### marker 1")
    except Exception as e:
        print("#####" + e)
    
    ## Add analysis info to data frame
    try:
        analysis.fdroop(spec_dict["substitutions"], spec_dict, df)
    except Exception as e:
        print(f'### fdroop Exception: {e}')
    # Save 
    with open(json_path, 'w') as json_file:
        json.dump(spec_dict,json_file, indent=2)
    
    if delete_src_data:
        os.remove(src_data_path)
        os.remove(src_json_path)
        
    return


def process_results_single_thread(
    spec_path: os.PathLike,
    temp_results_dir: os.PathLike,
    results_dir: os.PathLike,
    secs_to_remove_from_traces: float = 0,
    delete_src_data: bool = False,
    data_source_extension: str = ".psout",
):
    
    plotter = JRWFPlotter()
    analysis = Analysis()
    pdf_writer = PdfWriter()
    
    spec = None
    if not spec_path is None:
        while not os.path.exists(spec_path):
            time.sleep(30)
            logger.info(f"{spec_path} not found. Waiting 30 seconds...")
            
        logger.info(f"{spec_path} found.")
        spec = load_specs_from_csv(spec_path)

    # Results directory
    if spec_path is None:
        results_dir = os.path.join(results_dir,f"{get_date_time_str()}_results")
    # print(f"\nResults Directory located at: \n   {results_dir}")
    os.makedirs(results_dir,exist_ok=True)
    external_pdf_path = os.path.join(results_dir,"plot_summary.pdf")

    # Intermetiate temp results directory
    # print(f"\nTemp Results Directory located at: \n   {temp_results_dir}")
    os.makedirs(temp_results_dir,exist_ok=True)
    
    # List of sudies if spec is provided 
    study_list = []
    if spec_path:
        study_list = list(zip(list(spec["DIR"].values), list(spec["File_Name"].values)))
        num_studies = len(study_list)

    # While searching for results to process
    searching = True
    
    while searching:
        
        logger.info("Scanning for new results...")
    
        args = []
            
        for root, dirs, files in os.walk(temp_results_dir):
            
            for file in files:
                
                file_ext = os.path.splitext(file)[1]
                file_base_name = os.path.splitext(file)[0]
                
                
                if file_ext == data_source_extension:
                    
                    src_data_path = os.path.join(root, file_base_name + data_source_extension)
                    src_json_path = os.path.join(root, file_base_name + ".json")
                    
                    relative_path = os.path.relpath(root, temp_results_dir)
                    dst_results_dir = os.path.join(results_dir, relative_path)
                    os.makedirs(dst_results_dir, exist_ok=True)
                    
                    pkl_path = os.path.join(dst_results_dir, file_base_name + ".pkl")
                    png_path = os.path.join(dst_results_dir, file_base_name + ".png")
                    pdf_path = os.path.join(dst_results_dir, file_base_name + ".pdf")
                    json_path = os.path.join(dst_results_dir, file_base_name + ".json")
                    
                    # Arguments for processing study
                    process_results(
                        src_data_path,
                        src_json_path,
                        pkl_path,
                        png_path,
                        pdf_path,
                        json_path,
                        secs_to_remove_from_traces,
                        plotter,
                        analysis,
                        delete_src_data,
                    )
                    print("### marker 2")
                    # Appends to external results pdf
                    internal_pdf_file = open(pdf_path, 'rb')
                    print("### marker 22")
                    pdf_writer.append(fileobj=internal_pdf_file, pages=(0,1))
                    print("### marker 3")
                    # Remove study from study list if spec provided
                    if not spec is None:
                        while (relative_path,file_base_name) in study_list: 
                            study_list.remove((relative_path,file_base_name))
        print("### marker 4") 
        # Stop searching after one iteration if no spec provided
        if spec_path is None:
            logger.info("Stop Searching")
            searching = False
            
        # Stop searching after all studies removed from study list
        elif not study_list:
            searching = False
            logger.info("Processed all results")
        
        # wait between searches 
        time.sleep(30)
        
    # Save and close the pdf writer
    external_pdf = open(external_pdf_path, "wb")
    pdf_writer.write(external_pdf)
    pdf_writer.close()
    external_pdf.close()
    
    ## Plot the fdroop results
    fdroop_pdf_path = os.path.join(results_dir, "fdroop.pdf")
    
    try:
        analysis.plot_fdroop(fdroop_pdf_path)
    except Exception as e:
        print(f"### fdroop Exceptin: {e}")
            
    logger.warning("Exiting Process Results.")



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', "--output-spec-path", type=str, default=None)
    parser.add_argument('-t', "--temp-results-dir", type=str, default=None)
    parser.add_argument('-r', "--results-dir", type=str, default=None)
    parser.add_argument('-p', "--processes", type=int, default=1, help="(number of processes for multiprocessing)")
    parser.add_argument('-k', "--secs-to-remove", type=float, default=0, help="(seconds removed from begining of stored pkl)")
    parser.add_argument('-d', "--delete-src-data", type=bool, default=False, help="(1 = delete source data, 0 = keep source data)")
    parser.add_argument('-e', "--data-source", type=str, default=".psout",help="(.out or .pkl)")
    args = parser.parse_args()

    logger.info("Started: process_results.py")
    
    if args.temp_results_dir is None:
        print("Please Select Temp Directory:")
        temp_results_dir = prompt_for_directory_path(prompt_title="Select Temp Dir.", initial_dir="D:\\") #os.getcwd())
    else:
        temp_results_dir = args.temp_results_dir
        
    if args.results_dir is None:
        print("Please Select Results Directory:")
        results_dir = prompt_for_directory_path(prompt_title="Select Results Dir.", initial_dir="G:\\Junction_Rivers\\JRWF_PSCAD_Models") #os.getcwd())
    else:
        results_dir = args.results_dir
        
    spec_path = args.output_spec_path
    processes = args.processes
    secs_to_remove_from_traces = args.secs_to_remove
    delete_src_data = args.delete_src_data
    data_source_extension = args.data_source
    
    plotter = JRWFPlotter()
    analysis = Analysis()
    
    spec = None
    if not spec_path is None:
        while not os.path.exists(spec_path):
            time.sleep(30)
            logger.info(f"{spec_path} not found. Waiting 30 seconds...")
            
        logger.info(f"{spec_path} found.")
        spec = load_specs_from_csv(spec_path)

    # Results directory
    if spec_path is None:
        results_dir = os.path.join(results_dir,f"{get_date_time_str()}_results")
    print(f"\nResults Directory located at: \n   {results_dir}")
    os.makedirs(results_dir,exist_ok=True)

    # Intermetiate temp results directory
    # print(f"\nTemp Results Directory located at: \n   {temp_results_dir}")
    os.makedirs(temp_results_dir,exist_ok=True)
    
    # List of sudies if spec is provided 
    study_list = []
    if spec_path:
        study_list = list(zip(list(spec["DIR"].values), list(spec["File_Name"].values)))
        num_studies = len(study_list)

    # While searching for results to process
    searching = True
    while searching:
        
        logger.info("Scanning for new results...")
    
        args = []
            
        for root, dirs, files in os.walk(temp_results_dir):
            
            for file in files:
                
                file_ext = os.path.splitext(file)[1]
                file_base_name = os.path.splitext(file)[0]
                
                
                if file_ext == data_source_extension:
                    
                    src_data_path = os.path.join(root, file_base_name + data_source_extension)
                    src_json_path = os.path.join(root, file_base_name + ".json")
                    
                    relative_path = os.path.relpath(root, temp_results_dir)
                    dst_results_dir = os.path.join(results_dir, relative_path)
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
                        analysis,
                        delete_src_data,
                    ])
                    
                    # Remove study from study list if spec provided
                    if not spec is None:
                        while (relative_path,file_base_name) in study_list: 
                            study_list.remove((relative_path,file_base_name))
        
        if args: 
            remaining_str = f"({len(study_list)} of {num_studies} remaining)" if spec_path else ""
            logger.info(f"Processing {len(args)} Studies... {remaining_str}")
   
            # Multiprocess the processing of results
            with multiprocessing.Pool(processes=processes) as pool:
                pool.starmap(process_results, args)
        
        # Stop searching after one iteration if no spec provided
        if spec_path is None:
            logger.info("Stop Searching")
            searching = False
            
        # Stop searching after all studies removed from study list
        elif not study_list:
            searching = False
            logger.info("Processed all results")
        
        # wait between searches 
        time.sleep(30)
        
    logger.warning("Exiting Process Results.")