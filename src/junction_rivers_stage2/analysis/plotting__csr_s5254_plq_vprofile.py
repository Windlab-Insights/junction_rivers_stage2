import cmath
import math
import os
import re
from dataclasses import dataclass
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
import scipy.fft
import scipy.interpolate
from standard_mains import standard_plotting_main, find_filepath_beneath_dir
import json
import re
@dataclass
class FigureSet:
    scenarios: pd.DataFrame
    title: str
    filename: str

def get_file_path(directory, file_name):
    for root, _, files in os.walk(directory):
        if file_name in files:
            return os.path.join(root, file_name)
    return False

def init_matplotlib():
    plt.rcParams["axes.autolimit_mode"] = "round_numbers"
    plt.rcParams['axes.labelsize'] = 'medium'
    plt.rcParams['axes.titlesize'] = 'medium'
    plt.rcParams['axes.titleweight'] = 'bold'
    # plt.rcParams['axes.grid'] = True
    # plt.rcParams['axes.grid.axis'] = 'both'
    # plt.rcParams['axes.spines.right'] = False
    # plt.rcParams['axes.spines.top'] = False
    plt.rcParams['figure.titleweight'] = 'bold'
    plt.rcParams['figure.dpi'] = 300


def error_plot(tuning, figure_set: FigureSet, pkl_dir, output_dir, plot_length=7):
    df_subset = figure_set.scenarios
    fig_title = figure_set.title
    plot_filename = figure_set.filename


    init_matplotlib()
    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])


    fig, axs = plt.subplots(3, 1, gridspec_kw={'height_ratios': [2, 2, 1]})
    cm = 1/2.54
    # fig.set_size_inches(4*(1.3), 6*(1.3))
    fig.set_size_inches(10*cm, 18*cm)

    # Plot subtiltle.
    plt.suptitle(fig_title)
    plt.subplots_adjust(hspace=0.4, left=0.2, right=0.95, bottom=0.1)
    # plt.grid(color='lightgray', alpha=0.5, zorder=1)

    color = iter(plt.cm.inferno(np.linspace(0, 0.8, len(df_subset))))
    for _, scenario in df_subset.iterrows():

        dist_start = model_init_time + 1
        dist_duration = float(scenario['DistDuration'])
        dist_end = dist_start + dist_duration
        plot_start = dist_start - 0.5
        plot_end = dist_start + plot_length
        # plot_end = dist_start + 20
        print(f"Loading {scenario.File_Name}.pkl")
        pkl_path = get_file_path(pkl_dir,  scenario.File_Name + ".pkl")
        if os.path.isfile(pkl_path):
            c = next(color)
            data = pd.read_pickle(pkl_path)

            # Plot Ppoc Error.
            axs[0].title.set_text(r'$P_{\rm poc}$ Set-point Error')
            axs[0].grid(linestyle='--', color='lightgray', alpha=0.6)
            base_ppoc_mw = 400
            sig = (data['POC_P_MW'][plot_start:plot_end] - scenario['Init_Pwind_MW'] - scenario['Init_Pbess_MW'])
            sig = sig / base_ppoc_mw
            axs[0].plot(sig.index - dist_start, sig.values, linewidth=1, c=c)
            axs[0].set_ylabel(r'$P_{\rm poc}(t) - P_{\rm init} \quad {\rm (.pu)}$')
            axs[0].set_xlim(plot_start - dist_start, plot_end - dist_start)
            # axs[0].set_ylim(-10, 10)
            # axs[0].autoscale(enable=True, axis='x', tight=True)

            # Plot Qpoc Error.
            axs[1].title.set_text(r'$Q_{\rm poc}$ Set-point Error')
            base_qpoc_mvar = 158
            sig = (data['POC_Q_MVAr'][plot_start:plot_end] - data['Qref_droop_MVAr'][plot_start:plot_end])
            sig = sig / base_qpoc_mvar
            axs[1].grid(linestyle='--', color='lightgray', alpha=0.6)
            axs[1].plot(sig.index - dist_start, sig.values, linewidth=1, c=c)
            axs[1].set_ylabel(r'$Q_{\rm poc}(t) - Q_{\rm ref}(t) \quad {\rm (.pu)}$')
            axs[1].set_xlim(plot_start - dist_start, plot_end - dist_start)
            # # #     ax0.set_xscale('log')
            # axs[1].set_ylim(-10, 10)
            # axs[1].autoscale(enable=True, axis='x', tight=True)

            # Plot Trip Signals Error.
            axs[2].title.set_text(r'Trip Signals(solid) FRT Signals (dashed)')
            axs[2].grid(linestyle='--', color='lightgray', alpha=0.6)
            for sig_name in ['WT1_Trip', 'WT2_Trip', 'WT3_Trip', 'WT4_Trip']:
                sig = data[sig_name][plot_start:plot_end]
                axs[2].plot(sig.index - dist_start, sig.values, linewidth=1, c=c)
            for sig_name in ['WT1_FRT', 'WT2_FRT', 'WT3_FRT', 'WT4_FRT']:
                sig = data[sig_name][plot_start:plot_end]
                axs[2].plot(sig.index - dist_start, sig.values, '--', linewidth=1, c=c)

            # axs[2].set_ylabel(r'${\rm }$')
            axs[2].set_xlim(plot_start - dist_start, plot_end - dist_start)
            axs[2].set_ylim(-2, 2)
            axs[2].set_xlabel(r'Time from disturbance start (s)')
            # axs[1].autoscale(enable=True, axis='x', tight=True)
    outfile = os.path.join(output_dir, plot_filename)
    fig.align_labels()
    print("saving plot...")
    plt.savefig(outfile)
    plt.close()


def plotting_function(project_tuning, analysis_df, data_root_dir, output_dirpath):

    filename_prefix = "s5254_plq_vprofile_"

    min_flt_level = 523
    max_flt_level = 3138


    # FigureSet 1.
    filter = (analysis_df['Base_Title'].str.contains(r'^plq_id0[1-8]', regex=True))
    title = 'Results for PLQ Tests id: 1,2,4,5,6,7,8'
    filename = 'PLQ_Tests_ID_1_thru_8.png'
    fset1 = FigureSet(analysis_df[filter], title, filename_prefix + filename)

    # FigureSet 2.
    filter = (analysis_df['Base_Title'].str.contains(r'^plq_id(?:09|10|12)', regex=True))
    # filter = (analysis_df['Base_Title'].str.contains(r'^plq_id12', regex=True))
    title = 'PLQ Tests id: 9,10,12'
    filename = 'PLQ_Tests_ID_9_thru_12.png'
    fset2 = FigureSet(analysis_df[filter], title, filename_prefix + filename)

    filter = (analysis_df['Base_Title'].str.contains(r'^plq_id(?:0[1-9]|10)', regex=True))
    title = 'PLQ Tests id: 1,2,4,5,6,7,8,9,10'
    filename = 'PLQ_Tests_ID_1_thru_10.png'
    fset3 = FigureSet(analysis_df[filter], title, filename_prefix + filename)


    filter = (analysis_df['Base_Title'].str.contains(r'^plq_id(?:12)', regex=True))
    title = 'Results for PLQ Tests id: 12'
    filename = 'PLQ_Tests_ID_12.png'
    fset4 = FigureSet(analysis_df[filter], title, filename_prefix + filename)


    title = 'Results for PLQ CUO Tests'
    filename = 'PLQ_Tests_all.png'
    fset5 = FigureSet(analysis_df.head(12), title, filename_prefix + filename)

    error_plot(project_tuning, fset1, data_root_dir, output_dirpath, plot_length=9)
    error_plot(project_tuning, fset2, data_root_dir, output_dirpath, plot_length=9)
    error_plot(project_tuning, fset3, data_root_dir, output_dirpath, plot_length=9)
    error_plot(project_tuning, fset4, data_root_dir, output_dirpath, plot_length=15)
    error_plot(project_tuning, fset5, data_root_dir, output_dirpath, plot_length=9)





if __name__ == '__main__':
    standard_plotting_main(plotting_function=plotting_function)

