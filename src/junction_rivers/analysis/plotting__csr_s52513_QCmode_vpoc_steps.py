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


def qref_ppoc_error_plot(tuning, figure_set: FigureSet, pkl_dir, output_dir, plot_rise_time=True):
    df_subset = figure_set.scenarios
    fig_title = figure_set.title
    plot_filename = figure_set.filename
    plot_sig_label = r'$P_{\rm poc}(t) - P_{\rm init} \quad {\rm (.pu)}$'

    pbase_mw = float(tuning['SYS_Pbase_MW'])

    init_matplotlib()
    if df_subset['Ppoc_AEMO_Settling_Time'].empty:
        print("ppoc empty")
        return
    max_settle_time = max(df_subset[f'Ppoc_AEMO_Settling_Time'])
    print(f"MAX SETTLING Time {max_settle_time}")
    max_sig_settling_time = max(max(df_subset[f'Ppoc_AEMO_Settling_Time']), 1)
    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])

    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]})
    cm = 1 / 2.54
    # fig.set_size_inches(4*(1.3), 6*(1.3))
    fig.set_size_inches(12 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.17, right=0.95, bottom=0.1)

    ax = axs[0]
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    color = iter(plt.cm.inferno(np.linspace(0, 0.8, len(df_subset))))
    for _, scenario in df_subset.iterrows():

        time_steps = json.loads(scenario['Time_Steps'])
        vref_step_start = time_steps[0] + model_init_time
        vref_step_end = time_steps[1] + model_init_time

        plot_start = vref_step_start - (1 / 8) * max_sig_settling_time
        plot_end = plot_start + np.clip(4 * max_sig_settling_time, 1, 7.5)

        axs[0].title.set_text(fig_title)
        pkl_path = get_file_path(pkl_dir,  scenario.File_Name + ".pkl")
        if os.path.isfile(pkl_path):
            data = pd.read_pickle(pkl_path)

            pscad_signal = (
                        data['POC_P_MW'][plot_start:plot_end] - scenario['Init_Pwind_MW'] - scenario['Init_Pbess_MW'])
            pscad_signal = (pscad_signal / pbase_mw)
            # ax.cla()
            c = next(color)
            axs[0].plot(pscad_signal.index - vref_step_start, pscad_signal.values, linewidth=1, c=c)
        else:
            print(f"WARNING: No File found={pkl_path}")
            # raise Exception(f"No File found={pkl_path}")
    axs[0].set_ylabel(plot_sig_label)
    axs[0].autoscale(enable=True, axis='x', tight=True)
    axs[0].set_ylim(-0.05, 0.05)

    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    # ax.get_xaxis().set_visible(False)

    bdata = [list(df_subset[f'Ppoc_AEMO_Settling_Time'].values)]
    blabels = ['Settle-Time']

    axs[1].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[1].boxplot(bdata, vert=False, widths=0.5, manage_ticks=True)
    axs[1].set_yticklabels(blabels, rotation=45, fontsize=8)
    # axs[1].set_xlim(plot_start - vref_step_end, plot_end - vref_step_end)
    axs[1].set_xlim(axs[0].get_xlim())
    axs[1].set_xlabel('Time since Vpoc Disturbance (s)')

    # fig.align_labels()
    outfile = os.path.join(output_dir, plot_filename)
    plt.savefig(outfile)
    plt.close(fig)


def qref_vpoc_plot(tuning, figure_set: FigureSet, pkl_dir, output_dir, plot_rise_time=True):
    df_subset = figure_set.scenarios
    fig_title = figure_set.title
    plot_filename = figure_set.filename
    # plot_sig_label = r'$V_{\rm poc}(t) - V_{\rm ref,droop}(t) \quad {\rm (.pu)}$'
    plot_sig_label = r'$V_{\rm poc}(t)$'

    qbase_mvar = float(tuning['SYS_Qbase_MVAr'])
    droop_ratio = float(tuning['SYS_Vdroop_perc'])

    init_matplotlib()
    max_sig_settling_time = max(max(df_subset[f'Vpoc_AEMO_Settling_Time']), 1)
    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])

    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]})
    cm = 1 / 2.54
    # fig.set_size_inches(4*(1.3), 6*(1.3))
    fig.set_size_inches(12 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.17, right=0.95, bottom=0.1)

    ax = axs[0]
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    color = iter(plt.cm.inferno(np.linspace(0, 0.8, len(df_subset))))
    for _, scenario in df_subset.iterrows():

        time_steps = json.loads(scenario['Time_Steps'])
        vref_step_start = time_steps[0] + model_init_time
        vref_step_end = time_steps[1] + model_init_time

        plot_start = vref_step_start - (1 / 8) * max_sig_settling_time
        plot_end = plot_start + np.clip(4 * max_sig_settling_time, 1, 7.5)

        axs[0].title.set_text(fig_title)
        pkl_path = get_file_path(pkl_dir, scenario.File_Name + ".pkl")
        if os.path.isfile(pkl_path):
            data = pd.read_pickle(pkl_path)

            vpu = data['POC_V_pu'][plot_start:plot_end]
            # vref = data['Vref_pu'][plot_start:plot_end]
            qmvar = data['POC_Q_MVAr'][plot_start:plot_end]
            # qpu = qmvar / qbase_mvar
            # vref_droop_adjusted = vref - droop_ratio * qpu.clip(-1, 1)
            # pscad_signal = vpu - scenario['Vpoc_Settled_pu']
            pscad_signal = vpu
            # ax.cla()
            c = next(color)
            axs[0].plot(pscad_signal.index - vref_step_start, pscad_signal.values, linewidth=1, c=c)
        else:
            print(f"WARNING: No File found={pkl_path}")
            # raise Exception(f"No File found={pkl_path}")

    axs[0].set_ylabel(plot_sig_label)
    axs[0].autoscale(enable=True, axis='x', tight=True)
    # axs[0].set_ylim(-0.05, 0.05)

    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    # ax.get_xaxis().set_visible(False)

    bdata = [list(df_subset[f'Vpoc_AEMO_Settling_Time'].values)]
    blabels = ['Settle-Time']

    axs[1].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[1].boxplot(bdata, vert=False, widths=0.5, manage_ticks=True)
    axs[1].set_yticklabels(blabels, rotation=45, fontsize=8)
    # axs[1].set_xlim(plot_start - vref_step_end, plot_end - vref_step_end)
    axs[1].set_xlim(axs[0].get_xlim())
    axs[1].set_xlabel('Time since Vpoc Disturbance (s)')

    # fig.align_labels()
    outfile = os.path.join(output_dir, plot_filename)
    plt.savefig(outfile)
    plt.close(fig)


def qref_qpoc_error_plot(tuning, figure_set: FigureSet, pkl_dir, output_dir, plot_rise_time=False):
    df_subset = figure_set.scenarios
    fig_title = figure_set.title
    plot_filename = figure_set.filename
    plot_sig_label = r'$Q_{\rm poc}(t) - Q_{\rm ref}(t) \quad {\rm (.pu)}$'

    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])
    qbase_mvar = float(tuning['SYS_Qbase_MVAr'])
    droop_ratio = float(tuning['SYS_Vdroop_perc'])

    init_matplotlib()

    max_sig_settling_time = max(max(df_subset[f'Qpoc_AEMO_Settling_Time']), 1)

    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
    cm = 1 / 2.54
    # fig.set_size_inches(4*(1.3), 6*(1.3))
    fig.set_size_inches(12 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.17, right=0.95, bottom=0.1)

    ax = axs[0]
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)

    color = iter(plt.cm.inferno(np.linspace(0, 0.8, len(df_subset))))
    for _, scenario in df_subset.iterrows():
        # axs[0].title.text(fig_title)
        axs[0].title.set_text(fig_title)

        pkl_path = get_file_path(pkl_dir, scenario.File_Name + ".pkl")
        if os.path.isfile(pkl_path):
            data = pd.read_pickle(pkl_path)

            time_steps = json.loads(scenario['Time_Steps'])
            vref_step_start = time_steps[0] + model_init_time
            vref_step_end = time_steps[1] + model_init_time

            plot_start = vref_step_start - (1 / 8) * max_sig_settling_time
            plot_end = plot_start + np.clip(4 * max_sig_settling_time, 1, 10)

            #########################################################
            q_ctrl_mode = scenario['Control_Mode']
            if q_ctrl_mode == "VAR":
                qref_mvar = data["Qref_MVAr"][plot_start:plot_end]
            else:
                qref_mvar = data["Qref_droop_MVAr"][plot_start:plot_end]
            
            pscad_signal = (data['POC_Q_MVAr'][plot_start:plot_end] - qref_mvar)
            pscad_signal = pscad_signal / qbase_mvar
            # ax.cla()
            c = next(color)
            axs[0].plot(pscad_signal.index - vref_step_start, pscad_signal.values, linewidth=1, c=c)
        else:
            print(f"WARNING: No File found={pkl_path}")
            # raise Exception(f"No File found={pkl_path}")

    axs[0].set_ylabel(plot_sig_label)
    axs[0].autoscale(enable=True, axis='x', tight=True)
    # axs[0].set_ylim(-0.05, 0.05)

    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    # ax.get_xaxis().set_visible(False)

    if plot_rise_time:
        bdata = [list(df_subset[f'Qpoc_AEMO_Settling_Time'].values), list(df_subset[f'Qpoc_AEMO_Rise_Time'].values)]
        blabels = ['Settle-Time', 'Rise-Time']
    else:
        bdata = [list(df_subset[f'Qpoc_AEMO_Settling_Time'].values)]
        blabels = ['Settle-Time']


    axs[1].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[1].boxplot(bdata, vert=False, widths=0.5, manage_ticks=True)
    axs[1].set_yticklabels(blabels, rotation=45, fontsize=8)
    # axs[1].set_xlim(plot_start - vref_step_end, plot_end - vref_step_end)
    axs[1].set_xlim(axs[0].get_xlim())
    axs[1].set_xlabel('Time since Vpoc Disturbance (s)')

    # fig.align_labels()
    outfile = os.path.join(output_dir, plot_filename)
    plt.savefig(outfile)
    plt.close(fig)


def plotting_function(project_tuning, analysis_df, data_root_dir, output_dirpath):

    filename_prefix = "s52513_QCmode_vpoc_steps_"

    min_flt_level = 523
    max_flt_level = 3138
    FRT_MUST_ENGAGE_BY = 'Fault_Vpoc_below_80'

    # --- Ppoc Plots
    ppoc_settling_plots = []
    # exclude_categories = "|".join([re.escape(x) for x in ["Cont.Ctrl", "No-Windup[D.wc]", "No-Windup[U.wc]"]])
    # filter = (~analysis_df['Category'].str.contains(exclude_categories))
    filter = (analysis_df['Category'].isin(["d5%", "u5%"]))
    title = f'Q Control Mode Vpoc Step: Ppoc Error\n(All Tests)'
    filename = 'All_Ppoc_Error.png'
    ppoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    # # include_categories = "|".join([re.escape(x) for x in ["Sat(d5%)", "Sat(u5%)"]])
    # # filter = (analysis_df['Category'].str.contains(include_categories))
    # filter = (analysis_df['Category'].isin(["d5%", "u5%"]))
    # title = f'Q Control Mode Vpoc Step: Ppoc Error\n(5%-step No Saturation)'
    # filename = '5pc_NoSaturation_Ppoc_Error.png'
    # ppoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    # # include_categories = "|".join([re.escape(x) for x in ["PLQSat(d5%)", "PLQSat(u5%)"]])
    # # filter = (analysis_df['Category'].str.contains(include_categories))
    # filter = (analysis_df['Category'].isin(["Sat(d5%)", "Sat(u5%)"]))
    # title = f'Q Control Mode Vpoc Step: Ppoc Error\n(5%-step Saturation)'
    # filename = '5pc_Saturation_Ppoc_Error.png'
    # ppoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))


    for fset in ppoc_settling_plots:
        print(f"Saving {fset.title}")
        qref_ppoc_error_plot(project_tuning, fset, data_root_dir, output_dirpath)


    # --- Vpoc Plots
    vpoc_settling_plots = []
    # exclude_categories = "|".join([re.escape(x) for x in ["Cont.Ctrl", "No-Windup[D.wc]", "No-Windup[U.wc]"]])
    # filter = (~analysis_df['Category'].str.contains(exclude_categories))
    filter = (analysis_df['Category'].isin(["d5%", "u5%"]))
    title = 'Q Control Mode Q Control Mode: Vpoc\n(All Tests)'
    filename = 'All_Vpoc.png'
    vpoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    # # filter = (analysis_df['Category'].str.contains(include_categories))
    # filter = (analysis_df['Category'].isin(["d5%", "u5%"]))
    # title = 'Q Control Mode: Vpoc\n(5%-step No Saturation)'
    # filename = '5pc_No_Saturation_Vpoc.png'
    # vpoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    # filter = (analysis_df['Category'].isin(["Sat(d5%)", "Sat(u5%)"]))
    # title = f'Q Control Mode: Vpoc\n(5%-step Saturation)'
    # filename = '5pc_Saturation_Vpoc.png'
    # vpoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in vpoc_settling_plots:
        print(f"Saving {fset.title}")
        qref_vpoc_plot(project_tuning, fset, data_root_dir, output_dirpath, plot_rise_time=False)


    # --- Qpoc Plots
    qpoc_settling_plots = []
    filter = (analysis_df['Category'].isin(["d5%", "u5%"]))
    print(filter)
    title = 'Q Control Mode: Qpoc Error\n(All Tests)'
    filename = 'All_Qpoc_Error.png'
    qpoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    # # Naming of Sat here is confusing, it's based on mark/windlab definition. They now correspond to tests that
    # # don't saturate according to powerlink/TNSP definition.
    # # include_categories = "|".join([re.escape(x) for x in ["Sat(d5%)", "Sat(u5%)"]])
    # include_categories = ["d5%", "u5%"]
    # filter = (analysis_df['Category'].isin(include_categories))
    # title = 'Q Control Mode: Qpoc Error\n(5%-step No Saturation)'
    # filename = '5pc_No_Saturation_Qpoc_Error.png'
    # qpoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))


    # # include_categories = "|".join([re.escape(x) for x in ["PLQSat(d5%)", "PLQSat(u5%)"]])
    # include_categories = ["Sat(d5%)", "Sat(u5%)"]
    # filter = (analysis_df['Category'].isin(include_categories))
    # title = f'Q Control Mode: Qpoc Error\n(5% Step Saturation)'
    # filename = '5pc_Saturation_Qpoc_Error.png'
    # qpoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in qpoc_settling_plots:
        print(f"Saving {fset.title}")
        qref_qpoc_error_plot(project_tuning, fset, data_root_dir, output_dirpath)


if __name__ == '__main__':
    standard_plotting_main(plotting_function=plotting_function)

