from multiprocessing import Pool
import os
import subprocess
import pandas as pd

def save_filtered_df_as_csv(df: pd.DataFrame, key, col, path):


    filtered_df: pd.DataFrame = df.loc[df[col].str.contains(key)]

    filtered_df.to_csv(path)

    return path


analysis_scripts_dir = r"D:\grid_workspace\bungaban\analysis"
spec_dir = r"G:\Bungaban\Project_Parameters\Nordex_N163\specs_5000MVA_maxFL"
model_dir = r"G:\Bungaban\PSCAD_Models\Bungaban_Nordex_PSCAD_Model_v0.5"
results_dir = r"G:\Bungaban\PSCAD_Models\Bungaban_Nordex_PSCAD_Model_v0.5\results\2024_05_10_1758_results\2024_05_10_1759_results"
analysis_results_dir = os.path.join(results_dir,r"analysis")

# s5251
analysis_s5251_flatruns = [[
    "py", rf"{analysis_scripts_dir}\analysis__csr_s5251_flatruns.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-s", rf"{spec_dir}\specs__csr_s5251_flatruns.csv",
    "-t", rf"s5251",
    "-d", rf"{results_dir}\CSR_S5251",
    "-o", rf"{analysis_results_dir}\s5251_flatruns\s5251_flatruns_analysis.csv",
]]
# s5254
analysis_s5254_plq_profile = [[
    "py", rf"{analysis_scripts_dir}\plotting__csr_s5254_plq_vprofile.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-a", rf"{spec_dir}\specs__csr_s5254_plq_voltage_profile.csv",
    "-t", rf"s5254",
    "-d", rf"{results_dir}\CSR_S5254",
    "-o", rf"{analysis_results_dir}\s5254_plq_voltage_pyramid",
]]
# s5255 URES
analysis_s5255_ures = [
    [
    "py", rf"{analysis_scripts_dir}\analysis__csr_s5255_low_voltage_faults.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-s", rf"{spec_dir}\specs__csr_s5255_faults_ures.csv,{spec_dir}\specs__csr_s5255_faults_ohms.csv",
    "-t", rf"s5255",
    "-d", rf"{results_dir}\CSR_S5255",
    "-o", rf"{analysis_results_dir}\s5255_URES_OHMS\s5255_ures_analysis.csv",
],
# [
#     "py", rf"{analysis_scripts_dir}\plotting__csr_s5255_low_voltage_faults.py",
#     "-c", rf"{model_dir}\bbwf_tuning.json",
#     "-a", rf"{analysis_results_dir}\s5255_URES_OHMS\s5255_ures_analysis.csv",
#     "-t", rf"s5255",
#     "-d", rf"{results_dir}\CSR_S5255",
#     "-o", rf"{analysis_results_dir}\s5255_URES_OHMS",
# ]
]
# s5255 TOV
analysis_s5255_tov = [
    [
    "py", rf"{analysis_scripts_dir}\analysis__csr_s5255_high_voltage_faults.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-s", rf"{spec_dir}\specs__csr_s5255_faults_tov.csv",
    "-t", rf"s5255",
    "-d", rf"{results_dir}\CSR_S5255",
    "-o", rf"{analysis_results_dir}\s5255_TOV\s5255_tov_analysis_2.csv",
],
[
    "py", rf"{analysis_scripts_dir}\plotting__csr_s5255_high_voltage_faults.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-a", rf"{analysis_results_dir}\s5255_TOV\s5255_tov_analysis_2.csv",
    "-t", rf"s5255",
    "-d", rf"{results_dir}\CSR_S5255",
    "-o", rf"{analysis_results_dir}\s5255_TOV",
]]    
# s5257 Paretial Load Rejection
analysis_s5257_plr = [[
    "py", rf"{analysis_scripts_dir}\plotting__csr_s5257_plr.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-a", rf"{spec_dir}\specs__csr_s5257_partial_load_rejection.csv",
    "-t", rf"s52513",
    "-d", rf"{results_dir}\CSR_S5257",
    "-o", rf"{analysis_results_dir}\s55257_PLR",
]]
# s52513 VPOC Steps
analysis_s52513_vpoc_steps = [
    [
    "py", rf"{analysis_scripts_dir}\analysis__csr_s52513_vpoc_steps.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-s", rf"{spec_dir}\specs__csr_s52513_vpoc_steps.csv",
    "-t", rf"s52513",
    "-d", rf"{results_dir}\CSR_S52513",
    "-o", rf"{analysis_results_dir}\s52513_VPOC_STEPS\s52513_vpoc_steps_analysis.csv",
],
[
    "py", rf"{analysis_scripts_dir}\plotting__csr_s52513_vpoc_steps.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-a", rf"{analysis_results_dir}\s52513_VPOC_STEPS\s52513_vpoc_steps_analysis.csv",
    "-t", rf"s52513",
    "-d", rf"{results_dir}\CSR_S52513",
    "-o", rf"{analysis_results_dir}\s52513_VPOC_STEPS",
]]
# s52513 VRef Steps
analysis_s52513_vref_steps = [
    [
    "py", rf"{analysis_scripts_dir}\analysis__csr_s52513_vref_steps.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-s", rf"{spec_dir}\specs__csr_s52513_vref_steps.csv",
    "-t", rf"s52513",
    "-d", rf"{results_dir}\CSR_S52513",
    "-o", rf"{analysis_results_dir}\s52513_VREF_STEPS\s52513_vref_steps_analysis.csv",
],
[
    "py", rf"{analysis_scripts_dir}\plotting__csr_s52513_vref_steps.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-a", rf"{analysis_results_dir}\s52513_VREF_STEPS\s52513_vref_steps_analysis.csv",
    "-t", rf"s52513",
    "-d", rf"{results_dir}\CSR_S52513",
    "-o", rf"{analysis_results_dir}\s52513_VREF_STEPS",
]]
# s52513 QC VPOC VrefSteps
analysis_s52513_vpoc_qc_steps = [[
    "py", rf"{analysis_scripts_dir}\analysis__csr_s52513_QC_mode_qref_steps.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-s", rf"{spec_dir}\specs__csr_s52513_QC_mode_qref_steps.csv",
    "-t", rf"s52513",
    "-d", rf"{results_dir}\CSR_S52513",
    "-o", rf"{analysis_results_dir}\s52513_QC_model\s52513_QC_mode_qref_steps_analysis.csv",
],
[
    "py", rf"{analysis_scripts_dir}\plotting__csr_s52513_QCmode_qpoc_steps.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-a", rf"{analysis_results_dir}\s52513_QC_model\s52513_vpoc_steps_analysis.csv",
    "-t", rf"s52513",
    "-d", rf"{results_dir}\CSR_S52513",
    "-o", rf"{analysis_results_dir}\s52513_QC_model",
],
[
    "py", rf"{analysis_scripts_dir}\plotting__csr_s52513_QCmode_qref_steps.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-a", rf"{analysis_results_dir}\s52513_QC_model\s52513_QC_mode_qref_steps_analysis.csv",
    "-t", rf"s52513",
    "-d", rf"{results_dir}\CSR_S52513",
    "-o", rf"{analysis_results_dir}\s52513_QC_model",
]]
# s52513 OSC REJ
analysis_s52513_osc_rej = [[
    "py", rf"{analysis_scripts_dir}\plotting__csr_s52513_osc_rejection.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-a", rf"{spec_dir}\specs__wtg_dmat_bess_pzero_3215_osc_rejection.csv",
    "-t", rf"s52513",
    "-d", rf"{results_dir}\WTG_DMAT_BESS_PZERO_3215",
    "-o", rf"{analysis_results_dir}\s52513_OSC_REJ\WTG_DMAT_BESS_PZERO",
],
[
    "py", rf"{analysis_scripts_dir}\plotting__csr_s52513_osc_rejection.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-a", rf"{spec_dir}\specs__wtg_dmat_bess_pmax_3215_osc_rejection.csv",
    "-t", rf"s52513",
    "-d", rf"{results_dir}\WTG_DMAT_BESS_PMAX_3215",
    "-o", rf"{analysis_results_dir}\s52513_OSC_REJ\WTG_DMAT_BESS_PMAX",
],
[
    "py", rf"{analysis_scripts_dir}\plotting__csr_s52513_osc_rejection.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-a", rf"{spec_dir}\specs__wtg_dmat_bess_pmin_3215_osc_rejection.csv",
    "-t", rf"s52513",
    "-d", rf"{results_dir}\WTG_DMAT_BESS_PMIN_3215",
    "-o", rf"{analysis_results_dir}\s52513_OSC_REJ\WTG_DMAT_BESS_PMIN",
],
[
    "py", rf"{analysis_scripts_dir}\plotting__csr_s52513_osc_rejection.py",
    "-c", rf"{model_dir}\bbwf_tuning.json",
    "-a", rf"{spec_dir}\specs__bess_dmat_wtg_pzero_3215_osc_rejection.csv",
    "-t", rf"s52513",
    "-d", rf"{results_dir}\BESS_DMAT_WTG_PZERO_3215",
    "-o", rf"{analysis_results_dir}\s52513_OSC_REJ\BESS_DMAT_WTG_PZERO",
]]


analysis_layout = [
    # analysis_s5251_flatruns,
    # analysis_s5254_plq_profile,
    analysis_s5255_ures,
    # analysis_s5255_tov,
    # analysis_s5257_plr,
    # analysis_s52513_vpoc_qc_steps,
    # analysis_s52513_vpoc_steps,
    # analysis_s52513_vref_steps,
    # analysis_s52513_osc_rej,
]


def wrapper(args_collection):
    for args in args_collection:
        subprocess.run(args)

if __name__ == '__main__':

    with Pool(1) as pool:
        pool.map(wrapper, analysis_layout)

