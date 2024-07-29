import json
import os
import argparse
import pandas as pd
from rengen.utils.gui_utils import prompt_for_multiple_filepaths, prompt_for_directory_path, \
    prompt_for_saveas_filepath, prompt_for_open_filepath, save_csv_check_write_access
from rengen.spec.spec import load_specs_from_csv


def convert_to_absolute_path(path):
    if not os.path.isabs(path):
        # Path is relative, convert it to absolute
        absolute_path = os.path.abspath(path)
        return absolute_path
    else:
        # Path is already absolute
        return path


def find_filepath_beneath_dir(directory, file_name):
    for dir_path, _, files in os.walk(directory):
        if file_name in files:
            return os.path.join(dir_path, file_name)
    return None


def standard_analysis_main(file_name_prefix, per_scenario_analysis, summary_analysis=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--tuning-file', help='Paths to Project Tuning File', default=None)
    parser.add_argument('-s', '--spec-files', nargs="+", help='Paths to Flatrun Specs csvs', default=None)
    parser.add_argument('-t', '--test-prefix', nargs="+", help='Filename prefix for analysis.',
                        default=file_name_prefix)
    parser.add_argument('-d', '--data-root-dir', help='Parent folder of .pkl files.', default=None)
    parser.add_argument('-o', '--output', help='Output csv save file path', default=None)
    args = parser.parse_args()

    current_dir = os.getcwd()

    if args.tuning_file is None:
        print("Please Select Project Tuning File")
        tuning_file_path = prompt_for_open_filepath(prompt_title="Select Tuning File.")
        if tuning_file_path is None:
            exit()
    else:
        tuning_file_path = convert_to_absolute_path(args.tuning_file)

    with open(tuning_file_path) as f:
        project_tuning = json.load(f)

    if args.spec_files is None:
        print("Please Select Specs For Analysis.")
        spec_file_paths = prompt_for_multiple_filepaths(prompt_title="Select Spec Files.")
        if spec_file_paths is None:
            exit()
    else:
        print("args.spec_files: ",str(args.spec_files))
        spec_files_raw_paths = str(args.spec_files).strip("[] '").split(",")
        print("spec_files_raw_paths: ", spec_files_raw_paths)
        spec_file_paths = [convert_to_absolute_path(file) for file in spec_files_raw_paths]
        

    print("The following specs have been selected to Analysis:")
    for file_path in spec_file_paths:
        print(file_path)

    spec = load_specs_from_csv(spec_file_paths)

    print(f"Doing analysis over {args.test_prefix}")

    if args.data_root_dir is None:
        print("Please Select Data Root Directory for Analysis")
        data_root_dir = prompt_for_directory_path("Select Data(.pkl) Directory for Analysis", initial_dir=current_dir)
        if data_root_dir is None:
            exit()
    else:
        data_root_dir = convert_to_absolute_path(args.data_root_dir)

    print(f"Data Root Directory Selected: {data_root_dir}")

    if args.output is None:
        print("Please Select File-Path for Output CSV File.")
        output_filepath = prompt_for_saveas_filepath("Location to Save Output .csv File",
                                              initial_dir=current_dir,
                                              initialfile="per_scenario_analysis", defaultextension=".csv"
                                              )
        if output_filepath is None:
            exit()
    else:
        output_filepath = os.path.abspath(args.output)
        os.makedirs(os.path.dirname(output_filepath),exist_ok=True)

        

    # os.makedirs(os.path.dirname(output_filepath),exist_ok=True)

    print(f"Output Save File path Selected: {output_filepath}\n\tdir: {os.path.dirname(output_filepath)}")



    # Starting Analysis....
    # List of strings to match the start of 'File_Name'
    if isinstance(args.test_prefix, str):
        prefix_list = [args.test_prefix]

    else:
        prefix_list = args.test_prefix

    if args.test_prefix == ["None"]:
        prefix_list = None

    print(args.test_prefix)

    # Filter the DataFrame
    columns_to_drop = ['_2WorkerPartition', '_3WorkerPartition', '_4WorkerPartition', '_5WorkerPartition']
    columns_to_drop_existing = [col for col in columns_to_drop if col in spec.columns]
    if prefix_list:
        filtered_df = spec[spec['File_Name'].str.startswith(tuple(prefix_list))].drop(columns=columns_to_drop_existing)
    else:
        filtered_df = spec.drop(columns=columns_to_drop_existing)

    # Output Dataframe.
    per_scenario_analysis_df = pd.DataFrame.copy(filtered_df, deep=True).dropna(axis=1, how='all')

    # Do Some Analysis.
    for idx, scenario in enumerate(per_scenario_analysis_df.iloc):

        # Find Files Under Results Directory.
        scenario_pkl_filepath = find_filepath_beneath_dir(data_root_dir, f"{scenario['File_Name']}.pkl")
        scenario_metadata_filepath = find_filepath_beneath_dir(data_root_dir, f"{scenario['File_Name']}.json")
        if None not in [scenario_pkl_filepath, scenario_metadata_filepath]:
            print(f"Data exists for {scenario['File_Name']}.")

            data = pd.read_pickle(scenario_pkl_filepath)
            analysis_output = per_scenario_analysis(project_tuning, scenario, data)

            for k, v in analysis_output.items():
                filtered_df.loc[idx, k] = v


            # Check that scenario matches metadata.
            # metadata = pd.read_json(scenario_metadata_filepath).dropna(axis=1, how='all')
            #
            # for meta_key, meta_values in metadata.items():
            #     if scenario[meta_key] != meta_values:
            #         print(meta_key, meta_values)
            #         print(f"No Meta-data match for spec: {scenario['File_Name']}.")
            #     else:
            #         print(f"Meta-data matched: {scenario['File_Name']}.")

        else:
            print(f"No Data for {scenario['File_Name']}.")
    
    print("Saving df to csv... ")
    save_csv_check_write_access(filtered_df, output_filepath)

    if summary_analysis:
        summary_maps = summary_analysis(project_tuning, filtered_df, data_root_dir)
        # use the first rows keys as the columns for the whole dataframe.
        summary_df = pd.DataFrame(summary_maps, columns=list(summary_maps[0].keys()))
        save_csv_check_write_access(summary_df, f"{os.path.splitext(output_filepath)[0]}_summary.csv")



def standard_plotting_main(plotting_function):
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--tuning-file', help='Paths to Project Tuning File', default=None)
    parser.add_argument('-a', '--analysis-file', help='Paths to Analysis csv', default=None)
    parser.add_argument('-d', '--data-root-dir', help='Parent folder of .pkl files.', default=None)
    parser.add_argument('-t', '--test-prefix', nargs="+", help='Filename prefix for plotting.',
                        default=None)
    parser.add_argument('-o', '--output-dir', help='Plot output directory', default=None)
    args = parser.parse_args()

    current_dir = os.getcwd()

    if args.tuning_file is None:
        print("Please Select Project Tuning File")
        tuning_file_path = prompt_for_open_filepath(prompt_title="Select Tuning File.")
        if tuning_file_path is None:
            exit()
    else:
        tuning_file_path = convert_to_absolute_path(args.tuning_file)

    with open(tuning_file_path) as f:
        project_tuning = json.load(f)

    if args.analysis_file is None:
        print("Please Select Analysis File needed for plotting.")
        analysis_file = prompt_for_open_filepath(prompt_title="Select Analysis File")
        if analysis_file is None:
            exit()
    else:
        analysis_file = convert_to_absolute_path(args.analysis_file)

    print(f"The following Analysis file has been selected: {analysis_file}")

    spec = pd.read_csv(analysis_file)
    # print("\nspec",spec.head())

    

    if args.data_root_dir is None:
        print("Please Select Data Root Directory for Analysis")
        data_root_dir = prompt_for_directory_path("Select Data(.pkl) Directory for Analysis", initial_dir=current_dir)
        if data_root_dir is None:
            exit()
    else:
        data_root_dir = convert_to_absolute_path(args.data_root_dir)

    print(f"Data Root Directory Selected: {data_root_dir}")

    if args.output_dir is None:
        print("Please Select Directory-Path for Output Plots.")
        output_dirpath = prompt_for_directory_path("Directory to save plotting output.", initial_dir=current_dir)
        if output_dirpath is None:
            exit()
    else:
        output_dirpath = os.path.abspath(args.output_dir)


    if not os.path.exists(output_dirpath):
        os.makedirs(output_dirpath, exist_ok=True)


    print(f"Plot output directory path selected: {output_dirpath}")

    # Starting Plotting....
    # List of strings to match the start of 'File_Name'
    if args.test_prefix:
        if isinstance(args.test_prefix, str):
            prefix_list = [args.test_prefix]
        else:
            prefix_list = args.test_prefix
        
        
        filtered_df = spec[spec['File_Name'].str.startswith(str(prefix_list[-1]))]
    else:
        filtered_df = spec

    # Filter the DataFrame
    columns_to_drop = ['_2WorkerPartition', '_3WorkerPartition', '_4WorkerPartition', '_5WorkerPartition']
    columns_to_drop_existing = [col for col in columns_to_drop if col in spec.columns]
    filtered_df = filtered_df.drop(columns=columns_to_drop_existing)



    plotting_function(project_tuning, spec, data_root_dir, output_dirpath)
