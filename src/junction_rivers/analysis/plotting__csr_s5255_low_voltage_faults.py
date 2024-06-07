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

class PscadFaultType:
    NO_FAULT = 0
    A_TO_G = 1
    B_TO_G = 2
    C_TO_G = 3
    AB_TO_G = 4
    AC_TO_G = 5
    BC_TO_G = 6
    ABC_TO_G = 7
    A_TO_B = 8
    A_TO_C = 9
    B_TO_C = 10
    A_TO_B_TO_C = 11


@dataclass
class FigureSet:
    scenarios: pd.DataFrame
    title: str
    filename: str

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


def healthy_phase_box_and_whisker(tuning, figure_set: FigureSet, pkl_dir, output_dir):
    df_subset = figure_set.scenarios
    fig_title = figure_set.title
    plot_filename = figure_set.filename
    plot_sig_label = r'Max. In-fault POC Voltage (.pu)'

    min_flt_level = 523
    max_flt_level = 3138

    init_matplotlib()

    filter1ph_min = (df_subset['Fault_Type'] == 1) & (df_subset['Init_Fault_MVA'] == min_flt_level)
    filter1ph_max = (df_subset['Fault_Type'] == 1) & (df_subset['Init_Fault_MVA'] == max_flt_level)
    filter2ph_min = (df_subset['Fault_Type'] == 4) & (df_subset['Init_Fault_MVA'] == min_flt_level)
    filter2ph_max = (df_subset['Fault_Type'] == 4) & (df_subset['Init_Fault_MVA'] == max_flt_level)
    filterltl_min = (df_subset['Fault_Type'] == 8) & (df_subset['Init_Fault_MVA'] == min_flt_level)
    filterltl_max = (df_subset['Fault_Type'] == 8) & (df_subset['Init_Fault_MVA'] == max_flt_level)
    filter3phg_min = (df_subset['Fault_Type'] == 7) & (df_subset['Init_Fault_MVA'] == min_flt_level)
    filter3phg_max = (df_subset['Fault_Type'] == 7) & (df_subset['Init_Fault_MVA'] == max_flt_level)
    
    # Healthy Phase pu
    bdata_min_fault_lvl_healthy_phase = [
        list(df_subset[filter1ph_min]['Fault_Max_Healthy_Phase_pu'].values),
        list(df_subset[filter2ph_min]['Fault_Max_Healthy_Phase_pu'].values),
        list(df_subset[filterltl_min]['Fault_Max_Healthy_Phase_pu'].values),
        # list(df_subset[filter3phg_min]['Fault_Max_Healthy_Phase_pu'].values),  
    ]
    bdata_max_fault_lvl_healthy_phase = [
        list(df_subset[filter1ph_max]['Fault_Max_Healthy_Phase_pu'].values),
        list(df_subset[filter2ph_max]['Fault_Max_Healthy_Phase_pu'].values),
        list(df_subset[filterltl_max]['Fault_Max_Healthy_Phase_pu'].values),
        # list(df_subset[filter3phg_max]['Fault_Max_Healthy_Phase_pu'].values),   
    ]
    blabels = [
        'SLG',
        '2LG',
        'LL',
        # '3PG',
    ]
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[0].boxplot(bdata_min_fault_lvl_healthy_phase, widths=0.5, manage_ticks=True)
    axs[0].set_xticklabels(blabels)
    axs[0].set_ylabel(plot_sig_label)
    axs[0].set_title('In Fault Healthy Phase Voltages \n (Minimum Fault Level [523MVA])')

    axs[1].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[1].boxplot(bdata_max_fault_lvl_healthy_phase, widths=0.5, manage_ticks=True)
    axs[1].set_xticklabels(blabels)
    axs[1].set_ylabel(plot_sig_label)
    axs[1].set_title('In Fault Healthy Phase Voltages \n (Maximum Fault Level [3138MVA])')
    
    plt.tight_layout()
    outfile = os.path.join(output_dir, plot_filename)
    plt.savefig(outfile)
    plt.close(fig)
    
    in_fault_V_maxF = df_subset[filterltl_max]['Fault_Max_Healthy_Phase_pu'].values
    print(f'Max In_Fault V = {np.max(in_fault_V_maxF)}')
    in_fault_V_minF = df_subset[filterltl_min]['Fault_Max_Healthy_Phase_pu'].values
    print(f'Min In_Fault V = {np.max(in_fault_V_minF)}')

def max_phase_voltages_post_fault(tuning, figure_set: FigureSet, pkl_dir, output_dir):
    df_subset = figure_set.scenarios
    fig_title = figure_set.title
    plot_filename = figure_set.filename
    plot_sig_label = r'Max. Post-fault POC Voltage (.pu)'
    
    min_flt_level = 523
    max_flt_level = 3138

    init_matplotlib()

    filter1ph_min = (df_subset['Fault_Type'] == 1) & (df_subset['Init_Fault_MVA'] == min_flt_level)
    filter1ph_max = (df_subset['Fault_Type'] == 1) & (df_subset['Init_Fault_MVA'] == max_flt_level)
    filter2ph_min = (df_subset['Fault_Type'] == 4) & (df_subset['Init_Fault_MVA'] == min_flt_level)
    filter2ph_max = (df_subset['Fault_Type'] == 4) & (df_subset['Init_Fault_MVA'] == max_flt_level)
    filterltl_min = (df_subset['Fault_Type'] == 8) & (df_subset['Init_Fault_MVA'] == min_flt_level)
    filterltl_max = (df_subset['Fault_Type'] == 8) & (df_subset['Init_Fault_MVA'] == max_flt_level)
    filter3phg_min = (df_subset['Fault_Type'] == 7) & (df_subset['Init_Fault_MVA'] == min_flt_level)
    filter3phg_max = (df_subset['Fault_Type'] == 7) & (df_subset['Init_Fault_MVA'] == max_flt_level)
    
    # Psst fault V
    bdata_min_fault_lvl_post_fault = [
        list(df_subset[filter1ph_min]['Post_Fault_Max_V_pu'].values),
        list(df_subset[filter2ph_min]['Post_Fault_Max_V_pu'].values),
        list(df_subset[filterltl_min]['Post_Fault_Max_V_pu'].values),
        # list(df_subset[filter3phg_min]['Post_Fault_Max_V_pu'].values),  
    ]
    bdata_max_fault_lvl_post_fault = [
        list(df_subset[filter1ph_max]['Post_Fault_Max_V_pu'].values),
        list(df_subset[filter2ph_max]['Post_Fault_Max_V_pu'].values),
        list(df_subset[filterltl_max]['Post_Fault_Max_V_pu'].values),
        # list(df_subset[filter3phg_max]['Post_Fault_Max_V_pu'].values),   
    ]
    blabels = [
        'SLG',
        '2LG',
        'LL',
        # '3PG',
    ]
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[0].boxplot(bdata_min_fault_lvl_post_fault, widths=0.5, manage_ticks=True)
    axs[0].set_xticklabels(blabels)
    axs[0].set_ylabel(plot_sig_label)
    axs[0].set_title('Post fault Phase Voltages \n (Minimum Fault Level [523MVA])')

    axs[1].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[1].boxplot(bdata_max_fault_lvl_post_fault, widths=0.5, manage_ticks=True)
    axs[1].set_xticklabels(blabels)
    axs[1].set_ylabel(plot_sig_label)
    axs[1].set_title('Post fault Phase Voltages \n (Maximum Fault Level [3138MVA])')
    
    plt.tight_layout()
    outfile = os.path.join(output_dir, plot_filename)
    plt.savefig(outfile)
    plt.close(fig)
    
    

def healthy_phase_vs_fault_depth(tuning, figure_set: FigureSet, pkl_dir, output_dir):
    df_subset = figure_set.scenarios
    fig_title = figure_set.title
    plot_filename = figure_set.filename
    plot_sig_label = r'Post-fault Max. POC Voltage (.pu)'

    init_matplotlib()

    # fig, ax = plt.subplots(1, 1, gridspec_kw={'height_ratios': [3, 1]})
    fig, ax = plt.subplots(1)
    ax.title.set_text(fig_title)
    ax.set_ylabel(plot_sig_label)
    ax.set_xlabel(r'During-Fault Min. Faulted Phase Voltage (.pu)')

    cm = 1 / 2.54
    # fig.set_size_inches(4*(1.3), 6*(1.3))
    fig.set_size_inches(12 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)

    filter1ph = (df_subset['Fault_Type'] == 1)
    filter2ph = (df_subset['Fault_Type'] == 4)
    filterltl = (df_subset['Fault_Type'] == 8)
    filter3phg = (df_subset['Fault_Type'] == 7)

    ax.scatter(df_subset[filter1ph]['Fault_Min_Faulted_Phase_pu'].values,
                df_subset[filter1ph]['Post_Fault_Max_V_pu'].values, label="SLG")
    ax.scatter(df_subset[filter2ph]['Fault_Min_Faulted_Phase_pu'].values,
                df_subset[filter2ph]['Post_Fault_Max_V_pu'].values, label="2LG")
    ax.scatter(df_subset[filter3phg]['Fault_Min_Faulted_Phase_pu'].values,
                df_subset[filter3phg]['Post_Fault_Max_V_pu'].values, label="3LG")
    ax.scatter(df_subset[filterltl]['Fault_Min_Faulted_Phase_pu'].values,
                df_subset[filterltl]['Post_Fault_Max_V_pu'].values, label="LL")
    ax.legend()
    outfile = os.path.join(output_dir, plot_filename)
    plt.savefig(outfile)
    plt.close(fig)


def settling_time_plot(tuning, figure_set: FigureSet, pkl_dir, output_dir, plot_rise_time=True):
    df_subset = figure_set.scenarios
    if len(df_subset) == 0:
        return

    init_matplotlib()

    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])

    fig_title = figure_set.title
    plot_filename = figure_set.filename

    pscad_sig_name = 'POC_Iq_pos_pu'
    plot_sig_label = r'$i_{q,{\rm poc}}(t) \quad {\rm (.pu)}$'

    step_start =  model_init_time + df_subset['Fault_Time']
    max_sig_settling_time = max(df_subset[f'Fault_iq_AEMO_Settling_Time'])

    if min(df_subset['Fault_Duration']) != max(df_subset['Fault_Duration']):
        raise ValueError("Capacitive Iq Injection Settling time plot requires figure set to pass in sets of the same Fault_Duration.")
    else:
        fault_duration = max(df_subset['Fault_Duration'])

    # plot_start = step_start
    # plot_end = step_start + (10 / 8) * max_sig_settling_time
    # plot_end = min(plot_end, step_start + fault_duration)
    # plot_end = max(plot_end, step_start + 0.1)

    plot_duration = np.clip((10 / 8) * max_sig_settling_time, 0.1, fault_duration)

    cm = 1 / 2.54
    # fig.set_size_inches(4*(1.3), 6*(1.3))
    # fig.set_size_inches(10 * cm, 8 * cm)
    if plot_rise_time:
        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]})

    cm = 1 / 2.54
    fig.set_size_inches(12 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)

    color = iter(plt.cm.inferno(np.linspace(0, 0.8, len(df_subset))))
    for _, scenario in df_subset.sort_values(by=['Fault_Vpoc_pu'], ascending=False).iterrows():
        axs[0].title.set_text(fig_title)

        plot_start = model_init_time + scenario.Fault_Time
        plot_end = plot_start + plot_duration
        pkl_path = find_filepath_beneath_dir(pkl_dir, scenario.File_Name + ".pkl")

        if os.path.isfile(pkl_path):
            data = pd.read_pickle(pkl_path)

            pscad_signal = data[pscad_sig_name][plot_start:plot_end]
            # axs[0].cla()
            c = next(color)
            axs[0].plot(pscad_signal.index - plot_start, pscad_signal.values, linewidth=1, c=c)
        else:
            print(f"{scenario.File_Name}.pkl not found under {pkl_dir}")

    axs[0].set_ylabel(plot_sig_label)
    # Remove xlabels as we plan on covering plot. Don't want labels to poke out the sides.
    axs[0].set_xticklabels([])
    axs[0].autoscale(enable=True, axis='x', tight=True)
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)

    if plot_rise_time:
        bdata = [list(df_subset[f'Fault_iq_AEMO_Settling_Time'].values),
                 list(df_subset[f'Fault_iq_AEMO_Rise_Time'].values),
                #  list(df_subset[f'Iq_Above_MAS_Time'].values),
                 ]
        # blabels = ['Settling-time', 'Rise-time', 'Above-MAS']
        blabels = ['Settling-time', 'Rise-time']
    else:
        bdata = [list(df_subset[f'Fault_iq_AEMO_Settling_Time'].values),
                #  list(df_subset[f'Iq_Above_MAS_Time'].values),
                 ]
        # blabels = ['Settling-time', 'Above-MAS']
        blabels = ['Settling-time']
    
    axs[1].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[1].boxplot(bdata, vert=False, widths=0.5, manage_ticks=True)
    axs[1].set_yticklabels(blabels, rotation=45, fontsize=8)
    axs[1].set_xlim(0, plot_duration)
    axs[1].set_xlim(axs[0].get_xlim())
    axs[1].set_xlabel('Time since fault start (s)')
    # fig.align_labels()
    outfile = os.path.join(output_dir, plot_filename)
    plt.savefig(outfile)
    plt.close(fig)
    # fig.show()





def lv_active_power_recovery_plot(tuning, figure_set: FigureSet, pkl_dir, output_dir, plot_rise_time=True):
    raw_df_subset = figure_set.scenarios
    df_subset = raw_df_subset[raw_df_subset["Init_Pwind_MW"] != 0]
    if len(df_subset) == 0:
        return

    init_matplotlib()

    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])

    fig_title = figure_set.title
    plot_filename = figure_set.filename

    # WTG IQ Settling Time and rise time plots

    plot_sig_label = "WTG Active Power Error (MW)"

    
    max_sig_settling_time = max(df_subset['WTGS_95pc_Recovery_Time'])

    if min(df_subset['Fault_Duration']) != max(df_subset['Fault_Duration']):
        raise ValueError("active power recovery plot requires figure set to pass in sets of the same Fault_Duration.")
    else:
        fault_duration = max(df_subset['Fault_Duration'])

    plot_duration = np.clip((10 / 8) * max_sig_settling_time, 0.1, fault_duration)

    cm = 1 / 2.54
    if plot_rise_time:
        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]})

    cm = 1 / 2.54
    fig.set_size_inches(12 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)

    color = iter(plt.cm.inferno(np.linspace(0, 0.8, len(df_subset))))
    for _, scenario in df_subset.sort_values(by=['Fault_Vpoc_pu'], ascending=False).iterrows():
        axs[0].title.set_text("WTG "+fig_title)

        plot_start = model_init_time + scenario['Fault_Time'] + scenario['Fault_Duration']
        plot_end = plot_start + plot_duration
        pkl_path = find_filepath_beneath_dir(pkl_dir, scenario.File_Name + ".pkl")

        if scenario['Init_Pwind_MW'] == 0:
            continue

        if os.path.isfile(pkl_path):
            data = pd.read_pickle(pkl_path)

            pscad_signal = data['WT1_P_MW'][plot_start:plot_end] - data['Pref_WT_MW'][plot_start:plot_end]
            pscad_signal = data['WT2_P_MW'][plot_start:plot_end] - data['Pref_WT_MW'][plot_start:plot_end]
            pscad_signal = data['WT3_P_MW'][plot_start:plot_end] - data['Pref_WT_MW'][plot_start:plot_end]
            pscad_signal = data['WT4_P_MW'][plot_start:plot_end] - data['Pref_WT_MW'][plot_start:plot_end]
            # axs[0].cla()
            c = next(color)
            axs[0].plot(pscad_signal.index - plot_start, pscad_signal.values, linewidth=1, c=c)
        else:
            print(f"{scenario.File_Name}.pkl not found under {pkl_dir}")

    axs[0].set_ylabel(plot_sig_label)
    # Remove xlabels as we plan on covering plot. Don't want labels to poke out the sides.
    axs[0].set_xticklabels([])
    axs[0].autoscale(enable=True, axis='x', tight=True)
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)


    bdata = [
        list(df_subset[df_subset["WTGS_95pc_Recovery_Time"] != 0][f'WTGS_95pc_Recovery_Time'].values),
    ]
    # blabels = ['Settling-time', 'Rise-time', 'Above-MAS']
    blabels = ['95%\nRecovery']

    
    axs[1].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[1].boxplot(bdata, vert=False, widths=0.5, manage_ticks=True)
    axs[1].set_yticklabels(blabels, rotation=45, fontsize=8)
    axs[1].set_xlim(0, plot_duration)
    axs[1].set_xlim(axs[0].get_xlim())
    axs[1].set_xlabel('Time since fault start (s)')
    # fig.align_labels()
    outfile = os.path.join(output_dir, "WTG_"+plot_filename)
    plt.savefig(outfile)
    plt.close(fig)
    # fig.show()

    # BESS PLOTS

    plot_sig_label = "BESS Active Power Error (MW)"

    raw_df_subset = figure_set.scenarios
    df_subset = raw_df_subset[raw_df_subset["Init_Pbess_MW"] != 0]

    max_sig_settling_time = max(df_subset['BESS_95pc_Recovery_Time'])

    if min(df_subset['Fault_Duration']) != max(df_subset['Fault_Duration']):
        raise ValueError("active power recovery plot requires figure set to pass in sets of the same Fault_Duration.")
    else:
        fault_duration = max(df_subset['Fault_Duration'])

    plot_duration = np.clip((10 / 8) * max_sig_settling_time, 0.1, fault_duration)

    cm = 1 / 2.54
    if plot_rise_time:
        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]})

    cm = 1 / 2.54
    fig.set_size_inches(12 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)

    color = iter(plt.cm.inferno(np.linspace(0, 0.8, len(df_subset))))
    for _, scenario in df_subset.sort_values(by=['Fault_Vpoc_pu'], ascending=False).iterrows():
        axs[0].title.set_text("BESS "+fig_title)

        plot_start = model_init_time + scenario['Fault_Time'] + scenario['Fault_Duration']
        plot_end = plot_start + plot_duration
        pkl_path = find_filepath_beneath_dir(pkl_dir, scenario.File_Name + ".pkl")

        if scenario['Init_Pbess_MW'] == 0:
            continue

        if os.path.isfile(pkl_path):
            data = pd.read_pickle(pkl_path)

            data = pd.read_pickle(pkl_path)

            pscad_signal = data['BESS1_P_MW'][plot_start:plot_end] - data['Pref_BESS_ea_MW'][plot_start:plot_end]
            pscad_signal = data['BESS2_P_MW'][plot_start:plot_end] - data['Pref_BESS_ea_MW'][plot_start:plot_end]
            pscad_signal = data['BESS3_P_MW'][plot_start:plot_end] - data['Pref_BESS_ea_MW'][plot_start:plot_end]
            pscad_signal = data['BESS4_P_MW'][plot_start:plot_end] - data['Pref_BESS_ea_MW'][plot_start:plot_end]
            # axs[0].cla()
            c = next(color)
            axs[0].plot(pscad_signal.index - plot_start, pscad_signal.values, linewidth=1, c=c)
        else:
            print(f"{scenario.File_Name}.pkl not found under {pkl_dir}")

    axs[0].set_ylabel(plot_sig_label)
    # Remove xlabels as we plan on covering plot. Don't want labels to poke out the sides.
    axs[0].set_xticklabels([])
    axs[0].autoscale(enable=True, axis='x', tight=True)
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)


    bdata = [
        list(df_subset[df_subset["BESS_95pc_Recovery_Time"] != 0][f'BESS_95pc_Recovery_Time'].values),
    ]
    # blabels = ['Settling-time', 'Rise-time', 'Above-MAS']
    blabels = ['95%\nRecovery']

    
    axs[1].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[1].boxplot(bdata, vert=False, widths=0.5, manage_ticks=True)
    axs[1].set_yticklabels(blabels, rotation=45, fontsize=8)
    axs[1].set_xlim(0, plot_duration)
    axs[1].set_xlim(axs[0].get_xlim())
    axs[1].set_xlabel('Time since fault start (s)')
    # fig.align_labels()
    outfile = os.path.join(output_dir, "BESS_"+plot_filename)
    plt.savefig(outfile)
    plt.close(fig)
    # fig.show()



def terminal_settling_time_plot(tuning, figure_set: FigureSet, pkl_dir, output_dir, plot_rise_time=True):
    df_subset = figure_set.scenarios
    if len(df_subset) == 0:
        return

    init_matplotlib()

    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])

    fig_title = figure_set.title
    plot_filename = figure_set.filename

    # WTG IQ Settling Time and rise time plots

    plot_sig_label = r'$i_{q,{\rm poc}}(t) \quad {\rm (.pu)}$'

    step_start =  model_init_time + df_subset['Fault_Time']
    max_sig_settling_time = max(df_subset[f'Fault_iq_AEMO_Settling_Time'])

    if min(df_subset['Fault_Duration']) != max(df_subset['Fault_Duration']):
        raise ValueError("Capacitive Iq Injection Settling time plot requires figure set to pass in sets of the same Fault_Duration.")
    else:
        fault_duration = max(df_subset['Fault_Duration'])

    plot_duration = np.clip((10 / 8) * max_sig_settling_time, 0.1, fault_duration)

    cm = 1 / 2.54
    if plot_rise_time:
        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]})

    cm = 1 / 2.54
    fig.set_size_inches(12 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)

    color = iter(plt.cm.inferno(np.linspace(0, 0.8, len(df_subset))))
    for _, scenario in df_subset.sort_values(by=['Fault_Vpoc_pu'], ascending=False).iterrows():
        axs[0].title.set_text("WTG "+fig_title)

        plot_start = model_init_time + scenario.Fault_Time
        plot_end = plot_start + plot_duration
        pkl_path = find_filepath_beneath_dir(pkl_dir, scenario.File_Name + ".pkl")

        if os.path.isfile(pkl_path):
            data = pd.read_pickle(pkl_path)

            pscad_signal = data['WT1_Iq_pu'][plot_start:plot_end]
            pscad_signal = data['WT2_Iq_pu'][plot_start:plot_end]
            pscad_signal = data['WT3_Iq_pu'][plot_start:plot_end]
            pscad_signal = data['WT4_Iq_pu'][plot_start:plot_end]
            # axs[0].cla()
            c = next(color)
            axs[0].plot(pscad_signal.index - plot_start, pscad_signal.values, linewidth=1, c=c)
        else:
            print(f"{scenario.File_Name}.pkl not found under {pkl_dir}")

    axs[0].set_ylabel(plot_sig_label)
    # Remove xlabels as we plan on covering plot. Don't want labels to poke out the sides.
    axs[0].set_xticklabels([])
    axs[0].autoscale(enable=True, axis='x', tight=True)
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)

    if plot_rise_time:
        bdata = [list(df_subset[f'wtg_max_fault_iq_rise_time'].values),
                 list(df_subset[f'wtg_max_fault_iq_settle_time'].values),
                #  list(df_subset[f'Iq_Above_MAS_Time'].values),
                 ]
        # blabels = ['Settling-time', 'Rise-time', 'Above-MAS']
        blabels = ['Settling-time', 'Rise-time']
    else:
        bdata = [list(df_subset[f'Fault_iq_AEMO_Settling_Time'].values),
                #  list(df_subset[f'Iq_Above_MAS_Time'].values),
                 ]
        # blabels = ['Settling-time', 'Above-MAS']
        blabels = ['Settling-time']
    
    axs[1].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[1].boxplot(bdata, vert=False, widths=0.5, manage_ticks=True)
    axs[1].set_yticklabels(blabels, rotation=45, fontsize=8)
    axs[1].set_xlim(0, plot_duration)
    axs[1].set_xlim(axs[0].get_xlim())
    axs[1].set_xlabel('Time since fault start (s)')
    # fig.align_labels()
    outfile = os.path.join(output_dir, "WTG_"+plot_filename)
    plt.savefig(outfile)
    plt.close(fig)
    # fig.show()

    # BESS PLOTS

    cm = 1 / 2.54
    if plot_rise_time:
        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
    else:
        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]})

    cm = 1 / 2.54
    fig.set_size_inches(12 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)

    color = iter(plt.cm.inferno(np.linspace(0, 0.8, len(df_subset))))
    for _, scenario in df_subset.sort_values(by=['Fault_Vpoc_pu'], ascending=False).iterrows():
        axs[0].title.set_text("BESS " +fig_title)

        plot_start = model_init_time + scenario.Fault_Time
        plot_end = plot_start + plot_duration
        pkl_path = find_filepath_beneath_dir(pkl_dir, scenario.File_Name + ".pkl")

        if os.path.isfile(pkl_path):
            data = pd.read_pickle(pkl_path)

            pscad_signal = data['BESS1_Iq_pu'][plot_start:plot_end]
            pscad_signal = data['BESS2_Iq_pu'][plot_start:plot_end]
            pscad_signal = data['BESS3_Iq_pu'][plot_start:plot_end]
            pscad_signal = data['BESS4_Iq_pu'][plot_start:plot_end]
            # axs[0].cla()
            c = next(color)
            axs[0].plot(pscad_signal.index - plot_start, pscad_signal.values, linewidth=1, c=c)
        else:
            print(f"{scenario.File_Name}.pkl not found under {pkl_dir}")

    axs[0].set_ylabel(plot_sig_label)
    # Remove xlabels as we plan on covering plot. Don't want labels to poke out the sides.
    axs[0].set_xticklabels([])
    axs[0].autoscale(enable=True, axis='x', tight=True)
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)

    if plot_rise_time:
        bdata = [list(df_subset[f'BESS_max_fault_iq_rise_time'].values),
                 list(df_subset[f'BESS_max_fault_iq_settle_time'].values),
                #  list(df_subset[f'Iq_Above_MAS_Time'].values),
                 ]
        # blabels = ['Settling-time', 'Rise-time', 'Above-MAS']
        blabels = ['Settling-time', 'Rise-time']
    else:
        bdata = [list(df_subset[f'Fault_iq_AEMO_Settling_Time'].values),
                #  list(df_subset[f'Iq_Above_MAS_Time'].values),
                 ]
        # blabels = ['Settling-time', 'Above-MAS']
        blabels = ['Settling-time']
    
    axs[1].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[1].boxplot(bdata, vert=False, widths=0.5, manage_ticks=True)
    axs[1].set_yticklabels(blabels, rotation=45, fontsize=8)
    axs[1].set_xlim(0, plot_duration)
    axs[1].set_xlim(axs[0].get_xlim())
    axs[1].set_xlabel('Time since fault start (s)')
    # fig.align_labels()
    outfile = os.path.join(output_dir, "BESS_"+plot_filename)
    plt.savefig(outfile)
    plt.close(fig)
    # fig.show()



def ppoc_error_settling_time_plot(tuning, figure_set, pkl_dir, output_dir):
    df_subset = figure_set.scenarios
    if len(df_subset) == 0:
        return

    init_matplotlib()
    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])
    ppoc_mw_base = float(tuning['SYS_Pbase_MW'])

    fig_title = figure_set.title
    plot_filename = figure_set.filename
    plot_sig_label = r'$P_{\rm poc}(t) - P_{\rm init} \quad {\rm (.pu)}$'

    if min(df_subset['Fault_Duration']) != max(df_subset['Fault_Duration']):
        raise ValueError("URES PPOC Error Settling time plot requires figure set to pass in sets of the same Fault_Duration.")
    else:
        fault_duration = max(df_subset['Fault_Duration'])

    # fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
    fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [4, 1]})
    cm = 1 / 2.54
    # fig.set_size_inches(4*(1.3), 6*(1.3))
    fig.set_size_inches(13 * cm, 13 * cm)
    fig.subplots_adjust(hspace=0, left=0.18, right=0.95, bottom=0.1)

    max_sig_settling_time = max(df_subset[f'Ppoc_First_95pc_Recovery_Time'])

    ax = axs[0]
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)

    color = iter(plt.cm.inferno(np.linspace(0, 0.8, len(df_subset))))
    for _, scenario in df_subset.sort_values(by=['Fault_Vpoc_pu'], ascending=False).iterrows():
        step_start = model_init_time + scenario['Fault_Time']
        step_end = step_start + fault_duration

        plot_start = step_end - (1 / 8) * max_sig_settling_time
        plot_end = plot_start + max((10 / 8) * max_sig_settling_time, 0.1)
        # print(plot_start, plot_end)

        axs[0].title.set_text(fig_title)
        pkl_path = os.path.join(pkl_dir, scenario.File_Name + ".pkl")
        if os.path.isfile(pkl_path):
            data = pd.read_pickle(pkl_path)

            pscad_signal = (data['POC_P_MW'][plot_start:plot_end] - scenario['Init_Pwind_MW'] - scenario['Init_Pbess_MW'])
            pscad_signal = (pscad_signal / ppoc_mw_base)
            # ax.cla()
            c = next(color)
            axs[0].plot(pscad_signal.index - step_end, pscad_signal.values, linewidth=1, c=c)
        else:
            print(f"Not Found: {pkl_path}")
    axs[0].set_ylabel(plot_sig_label)
    axs[0].autoscale(enable=True, axis='x', tight=True)
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[0].set_xticklabels([]) # Remove ticklabels as we plan on them.
    curr_ylims = axs[0].get_ylim()

    x_vals = np.linspace(plot_start- step_end, plot_end - step_end, 100)
    settle_line = -0.05 * ((scenario['Init_Pwind_MW'] + scenario['Init_Pbess_MW']) / ppoc_mw_base)
    y_vals = [settle_line] * len(x_vals)
    fill_vals = [np.sign(settle_line)] * len(x_vals)
    axs[0].axhline(settle_line, color='black', linestyle='--', linewidth='0.5')
    axs[0].axhline(0, color='black', linestyle='-', linewidth='0.5')
    axs[0].fill_between(x_vals, y_vals, fill_vals, color='yellow', alpha=0.05)

    axs[0].set_ylim(curr_ylims)
    
    bdata = [df_subset[f'Ppoc_First_95pc_Recovery_Time'].values]
    blabels = ['First\nWithin\n95%']
    
    axs[1].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    axs[1].boxplot(bdata, vert=False, widths=0.5, manage_ticks=True)
    axs[1].set_yticklabels(blabels, rotation=45, fontsize=8)
    # axs[1].set_xlim(plot_start - vref_step_end, plot_end - vref_step_end)
    axs[1].set_xlim(axs[0].get_xlim())
    axs[1].set_xlabel('Time since fault cleared (s)')

    # fig.align_labels()
    outfile = os.path.join(output_dir, plot_filename)
    plt.savefig(outfile)
    plt.close(fig)
    # fig.show()



#
# def ppoc_error_settling_time_plot(tuning, figure_set, pkl_dir, output_dir):
#     df_subset = figure_set.scenarios
#     if len(df_subset) == 0:
#         return
#
#     init_matplotlib()
#     model_init_time = float(tuning['TIME_Full_Init_Time_sec'])
#     ppoc_mw_base = float(tuning['SYS_Pbase_MW'])
#
#     fig_title = figure_set.title
#     plot_filename = figure_set.filename
#     plot_sig_label = r'$P_{\rm poc}(t) - P_{\rm init} \quad {\rm (.pu)}$'
#
#     if min(df_subset['Fault_Duration']) != max(df_subset['Fault_Duration']):
#         raise ValueError("URES PPOC Error Settling time plot requires figure set to pass in sets of the same Fault_Duration.")
#     else:
#         fault_duration = max(df_subset['Fault_Duration'])
#
#     fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]})
#     cm = 1 / 2.54
#     # fig.set_size_inches(4*(1.3), 6*(1.3))
#     fig.set_size_inches(12 * cm, 12 * cm)
#     plt.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)
#
#     max_sig_settling_time = max(df_subset[f'Ppoc_Within20MW_Recovery_Time'])
#
#     axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
#
#     color = iter(plt.cm.inferno(np.linspace(0, 0.8, len(df_subset))))
#     for _, scenario in df_subset.sort_values(by=['Fault_Vpoc_pu'], ascending=False).iterrows():
#         step_start = model_init_time + scenario['Fault_Time']
#         step_end = step_start + fault_duration
#
#         plot_start = step_end - (1 / 8) * max_sig_settling_time
#         plot_end = plot_start + (10 / 8) * max_sig_settling_time
#
#         axs[0].title.set_text(fig_title)
#         pkl_path = os.path.join(pkl_dir, scenario.File_Name + ".pkl")
#         if os.path.isfile(pkl_path):
#             data = pd.read_pickle(pkl_path)
#
#             pscad_signal = (data['POC_P_MW'][plot_start:plot_end] - scenario['Init_Pwind_MW'] - scenario['Init_Pbess_MW'])
#             pscad_signal = (pscad_signal / ppoc_mw_base)
#             # ax.cla()
#             c = next(color)
#             axs[0].plot(pscad_signal.index - step_end, pscad_signal.values, linewidth=1, c=c)
#     axs[0].set_ylabel(plot_sig_label)
#     axs[0].autoscale(enable=True, axis='x', tight=True)
#     axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
#     axs[0].set_xticklabels([]) # Remove ticklabels as we plan on them.
#
#     new_metric = np.minimum(df_subset[f'Ppoc_Within10MW_Recovery_Time'].values, df_subset[f'Ppoc_Within20MW_Recovery_Time'].values)
#     bdata = [list(new_metric), list(df_subset[f'Ppoc_Within20MW_Recovery_Time'].values)]
#     blabels = ['95% or 10MW', '95%\nrecovery']
#
#     axs[1].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
#     axs[1].boxplot(bdata, vert=False, widths=0.5, manage_ticks=True)
#     axs[1].set_yticklabels(blabels, rotation=45, fontsize=8)
#     # axs[1].set_xlim(plot_start - vref_step_end, plot_end - vref_step_end)
#     axs[1].set_xlim(axs[0].get_xlim())
#     axs[1].set_xlabel('Time since fault cleared (s)')
#
#     # fig.align_labels()
#     outfile = os.path.join(output_dir, plot_filename)
#     plt.savefig(outfile)
#     plt.close(fig)
#     # fig.show()
#
    

 

def terminal_low_voltage_iq_injection_plot(tuning, figure_set, pkl_dir, output_dir):

    
    df_subset = figure_set.scenarios
    if len(df_subset) == 0:
        return
    
    voltage_signals = [
        'WT1_Fault_V_pu',
        'WT2_Fault_V_pu',
        'WT3_Fault_V_pu',
        'WT4_Fault_V_pu',
        'BESS1_Fault_V_pu',
        'BESS2_Fault_V_pu',
        'BESS3_Fault_V_pu',
        'BESS4_Fault_V_pu',
    ]

    for voltage_signal in voltage_signals:
        df_subset = df_subset[df_subset[voltage_signal] < 0.8]


    init_matplotlib()
    fig_title = figure_set.title
    plot_filename = figure_set.filename

    def lvrt_iq_guide(k, u, v):
        return np.clip(k*(u - v), 0, 1)

    v_array = np.linspace(0, 1, 100)
    wtg_nas_perc = 3
    bess_nas_perc = 0.8
    iq_mas_values = [lvrt_iq_guide(0, 0.8, v) for v in v_array]
    wtg_iq_nas_values = [lvrt_iq_guide(wtg_nas_perc, 0.8, v) for v in v_array]
    bess_iq_nas_values = [lvrt_iq_guide(bess_nas_perc, 0.8, v) for v in v_array]
    iq_ass_values = [lvrt_iq_guide(4, 0.85, v) for v in v_array]

    figs, ax = plt.subplots(1)
    cm = 1 / 2.54
    figs.set_size_inches(12 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)

    ax.plot(v_array, iq_mas_values, color='yellow', label='mas[0%/%]')
    ax.plot(v_array, wtg_iq_nas_values, color='gray', label=f'nas[{wtg_nas_perc}%/%]')
    ax.plot(v_array, iq_ass_values, color='green', label='aas[4%/%]')
    ax.fill_between(v_array, wtg_iq_nas_values, iq_mas_values, color='yellow', alpha=0.1)
    ax.fill_between(v_array, iq_ass_values, wtg_iq_nas_values, color='gray', alpha=0.1)
    ax.fill_between(v_array, [2] * len(v_array), iq_ass_values, color='green', alpha=0.1)
    ax.grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    ax.scatter(df_subset['WT1_Fault_V_pu'], df_subset['WT1_Fault_iq_pu'],  edgecolor="black", facecolor='blue', alpha=0.6, zorder=2)
    ax.scatter(df_subset['WT2_Fault_V_pu'], df_subset['WT2_Fault_iq_pu'],  edgecolor="black", facecolor='blue', alpha=0.6, zorder=2)
    ax.scatter(df_subset['WT3_Fault_V_pu'], df_subset['WT3_Fault_iq_pu'],  edgecolor="black", facecolor='blue', alpha=0.6, zorder=2)
    ax.scatter(df_subset['WT4_Fault_V_pu'], df_subset['WT4_Fault_iq_pu'],  edgecolor="black", facecolor='blue', alpha=0.6, zorder=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.1, 1.5)
    ax.set_title("WTG "+fig_title)
    ax.set_ylabel(f"$i_q$ Injection [pu mcc]")
    ax.set_xlabel("WTG Terminal Settled Fault Voltage [pu]")
    ax.legend()
    figs.savefig(os.path.join(output_dir, "wtg_"+plot_filename))

    figs, ax = plt.subplots(1)
    cm = 1 / 2.54
    figs.set_size_inches(12 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)

    ax.plot(v_array, iq_mas_values, color='yellow', label='mas[0%/%]')
    ax.plot(v_array, bess_iq_nas_values, color='gray', label=f'nas[{bess_nas_perc}%/%]')
    ax.plot(v_array, iq_ass_values, color='green', label='aas[4%/%]')
    ax.fill_between(v_array, bess_iq_nas_values, iq_mas_values, color='yellow', alpha=0.1)
    ax.fill_between(v_array, iq_ass_values, bess_iq_nas_values, color='gray', alpha=0.1)
    ax.fill_between(v_array, [2] * len(v_array), iq_ass_values, color='green', alpha=0.1)
    ax.grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    ax.scatter(df_subset['BESS1_Fault_V_pu'], df_subset['BESS1_Fault_iq_pu'], edgecolor="black", facecolor='green', alpha=0.6, zorder=2)
    ax.scatter(df_subset['BESS2_Fault_V_pu'], df_subset['BESS2_Fault_iq_pu'], edgecolor="black", facecolor='green', alpha=0.6, zorder=2)
    ax.scatter(df_subset['BESS3_Fault_V_pu'], df_subset['BESS3_Fault_iq_pu'], edgecolor="black", facecolor='green', alpha=0.6, zorder=2)
    ax.scatter(df_subset['BESS4_Fault_V_pu'], df_subset['BESS4_Fault_iq_pu'], edgecolor="black", facecolor='green', alpha=0.6, zorder=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.1, 1.5)
    ax.set_title("BESS "+fig_title)
    ax.set_ylabel(f"$i_q$ Injection [pu mcc]")
    ax.set_xlabel("BESS Terminal Settled Fault Voltage [pu]")
    ax.legend()
    figs.savefig(os.path.join(output_dir, "bess_"+plot_filename))




def low_voltage_iq_injection_plot(tuning, figure_set, pkl_dir, output_dir):
    df_subset = figure_set.scenarios
    if len(df_subset) == 0:
        return

    init_matplotlib()
    fig_title = figure_set.title
    plot_filename = figure_set.filename

    def lvrt_iq_guide(k, u, v):
        return np.clip(k*(u - v), 0, 1)

    v_array = np.linspace(0, 1, 100)
    nas_perc = 0.8
    iq_mas_values = [lvrt_iq_guide(0, 0.8, v) for v in v_array]
    iq_nas_values = [lvrt_iq_guide(nas_perc, 0.8, v) for v in v_array]
    iq_ass_values = [lvrt_iq_guide(4, 0.85, v) for v in v_array]

    figs, ax = plt.subplots(1)
    cm = 1 / 2.54
    figs.set_size_inches(12 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)

    ax.plot(v_array, iq_mas_values, color='yellow', label='mas[0%/%]')
    ax.plot(v_array, iq_nas_values, color='gray', label=f'nas[{nas_perc}%/%]')
    ax.plot(v_array, iq_ass_values, color='green', label='aas[4%/%]')
    ax.fill_between(v_array, iq_nas_values, iq_mas_values, color='yellow', alpha=0.1)
    ax.fill_between(v_array, iq_ass_values, iq_nas_values, color='gray', alpha=0.1)
    ax.fill_between(v_array, [2] * len(v_array), iq_ass_values, color='green', alpha=0.1)
    ax.grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    ax.scatter(df_subset['Fault_Vpoc_pu'], df_subset['Delta_iq_mcc_pu'], edgecolor='black', facecolor='white', zorder=2)
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.1, 1.5)
    ax.set_title(fig_title)
    ax.set_ylabel(f"$i_q$ Injection [pu mcc]")
    ax.set_xlabel("POC Settled Fault Voltage [pu]")
    ax.legend()
    figs.savefig(os.path.join(output_dir, plot_filename))

    # ax.plot(v_array, iq_mas_values, ls='-', color='yellow', label='mas[0%/%]')
    # ax.plot(v_array, iq_mas_old_values, ls='-', color='gray', label='nas[3%/%]')
    # ax.plot(v_array, iq_ass_values, ls='-', color='green', label='aas[6%/%]')
    #
    # ax.fill_between(v_array, iq_mas_old_values, iq_mas_values, color='yellow', alpha=0.1)
    # ax.fill_between(v_array, iq_ass_values, iq_mas_old_values, color='gray', alpha=0.1)
    # ax.fill_between(v_array, [-2] * len(v_array), iq_ass_values, color='green', alpha=0.1)
    # ax.grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)




def rise_and_settle_times_against_fault_delta_v(tuning, figure_set, pkl_dir, output_dir):
    df_subset = figure_set.scenarios
    if len(df_subset) == 0:
        return

    init_matplotlib()
    fig_title = figure_set.title
    plot_filename = figure_set.filename


    figs, axs = plt.subplots(nrows=1, ncols=2)
    cm = 1 / 2.54
    figs.set_size_inches(24 * cm, 12 * cm)
    plt.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)

    print(axs)

    ax = axs[0]
    ax.axhline(0.040, linestyle='--', label='AAS[40ms]')
    ax.grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    ax.scatter(df_subset['Fault_Vpoc_Range_pu'], df_subset['Fault_iq_AEMO_Rise_Time'], edgecolor='black', facecolor='white', zorder=2)
    # ax.set_xlim(0, 1)
    # ax.set_ylim(-0.1, 1.5)
    ax.set_title(f"{fig_title}\nIq Rise-Time Summary")
    ax.set_ylabel(f"$i_q$ rise-time [secs]")
    ax.set_xlabel("Vpoc deviation during fault [pu]")
    ax.legend()

    ax = axs[1]
    ax.axhline(0.070, linestyle='--', label='AAS[70ms]')
    ax.grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)
    ax.scatter(df_subset['Fault_Vpoc_Range_pu'], df_subset['Fault_iq_AEMO_Settling_Time'], edgecolor='black', facecolor='white', zorder=2)
    # ax.set_xlim(0, 1)
    # ax.set_ylim(-0.1, 1.5)
    ax.set_title(f"{fig_title}\nIq Settle-Time Summary")
    ax.set_ylabel(f"$i_q$ settle-time [secs]")
    ax.set_xlabel("Vpoc deviation during fault [pu]")
    ax.legend()

    figs.savefig(os.path.join(output_dir, plot_filename))

def plotting_function(project_tuning, analysis_df, data_root_dir, output_dirpath):

    filename_prefix = "s5255_low_voltage_"

    min_flt_level = 523
    max_flt_level = 3138
    FRT_MUST_ENGAGE_BY = 'Fault_Vpoc_below_80'

    # --- Low Voltage IQ Injection plots
    iq_injection_fsets = []
    
    std_filter_ures = ((analysis_df['Category'] != 'Fault_Ohms') &
                  (analysis_df['Fault_Duration'] == 0.43) &
                  (analysis_df['Fault_iq_mcc_pu'] < 1.0) &
                  (analysis_df['BESS_Irms_over_190']== False) &
                  (analysis_df[FRT_MUST_ENGAGE_BY] == True))
    
    filter = (std_filter_ures &
              (analysis_df['Init_Fault_MVA'] == max_flt_level) & (analysis_df['Fault_Type'] == PscadFaultType.ABC_TO_G))
    title = "POC $i_q$ Injection: Max Fault Level [3138MVA]\nBalanced Under-voltage Faults"
    filename = 'POC_Balanced_Iq_Injection_MaxFaultLevel.png'
    iq_injection_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = (std_filter_ures &
              (analysis_df['Init_Fault_MVA'] == min_flt_level) & (analysis_df['Fault_Type'] == PscadFaultType.ABC_TO_G))
    title = "POC $i_q$ Injection: Min. Fault Level [523MVA]\nBalanced Under-voltage Faults"
    filename = 'POC_Balanced_Iq_Injection_MinFaultLevel.png'
    iq_injection_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = (std_filter_ures &
              (analysis_df['Init_Fault_MVA'] == max_flt_level) & (analysis_df['Fault_Type'] != PscadFaultType.ABC_TO_G))
    title = "POC $i_q$ Injection: Max Fault Level [3138MVA]\nUnbalanced Under-voltage Faults"
    filename = 'POC_Unbalanced_Iq_Injection_MaxFaultLevel.png'
    iq_injection_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = (std_filter_ures &
              (analysis_df['Init_Fault_MVA'] == min_flt_level) & (analysis_df['Fault_Type'] != PscadFaultType.ABC_TO_G))
    title = "POC $i_q$ Injection: Min. Fault Level [523MVA]\nUnbalanced Under-voltage Faults"
    filename = 'POC_Unbalanced_Iq_Injection_MinFaultLevel.png'
    iq_injection_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    
    filter = (std_filter_ures)
    title = "POC $i_q$ Injection: \nAll Non-resistive Faults"
    filename = 'POC_All_Non-resis_Iq_Injection.png'
    iq_injection_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    
    std_filter_resistive = ((analysis_df['Category'] == 'Fault_Ohms') &
                  (analysis_df['Fault_Duration'] == 0.43) &
                  (analysis_df['Fault_iq_mcc_pu'] < 1.0) &
                  (analysis_df[FRT_MUST_ENGAGE_BY] == True))
    
    filter = (std_filter_resistive & (analysis_df['Init_Fault_MVA'] == max_flt_level))
    title = "POC $i_q$ Injection: Max Fault Level [3138MVA]\nResistive Faults"
    filename = 'POC_Resistive_Iq_Injection_MaxFaultLevel.png'
    iq_injection_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = (std_filter_resistive & (analysis_df['Init_Fault_MVA'] == min_flt_level))
    title = "POC $i_q$ Injection: Min. Fault Level [523MVA]\nResistive Faults"
    filename = 'POC_Resistive_Iq_Injection_MinFaultLevel.png'
    iq_injection_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    
    std_filter = ((analysis_df['Fault_Duration'] == 0.43) &
                  (analysis_df['Fault_iq_mcc_pu'] < 1.0) &
                  (analysis_df[FRT_MUST_ENGAGE_BY] == True))

    filter = std_filter
    title = "POC $i_q$ Injection\nAll Under-voltage Faults"
    filename = 'POC_ALL_Iq_Injection.png'
    iq_injection_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in iq_injection_fsets:
        print(f"Saving {fset.title}")
        low_voltage_iq_injection_plot(project_tuning, fset, data_root_dir, output_dirpath)

    terminal_iq_injetion_fsets = []

    filter = std_filter
    title = "$i_q$ Injection at Terminals\nAll Under-voltage Faults"
    filename = 'Terminal_ALL_Iq_Injection.png'
    terminal_iq_injetion_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = (std_filter_ures)
    title = "$i_q$ Injection at Terminals: \nAll Non-resistive Faults"
    filename = 'Terminal_All_Non-resis_Iq_Injection.png'
    terminal_iq_injetion_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in terminal_iq_injetion_fsets:
        print(f"Saving {fset.title}")
        terminal_low_voltage_iq_injection_plot(project_tuning, fset, data_root_dir, output_dirpath)
        
        
    terminal_iq_injetion_fsets = []

    filter = std_filter
    title = "Terminal Capacitive Current Injection: \nAll Under-voltage Faults"
    filename = 'Terminal_ALL_Iq_Injection.png'
    terminal_iq_injetion_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = (std_filter_ures)
    title = "Terminal Capacitive Current Injection: \nAll Non-resistive Faults"
    filename = 'Terminal_All_Non-resis_Iq_Injection.png'
    terminal_iq_injetion_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in terminal_iq_injetion_fsets:
        print(f"Saving {fset.title}")
        terminal_settling_time_plot(project_tuning, fset, data_root_dir, output_dirpath, plot_rise_time=True)


    lv_active_power_fsets = []

    filter = std_filter
    title = "Fault Active Power: \nAll Under-voltage Faults"
    filename = 'Terminal_ALL_active_power.png'
    lv_active_power_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = (std_filter_ures)
    title = "Fault Active Power: \nAll Non-resistive Faults"
    filename = 'Terminal_All_Non-resis_active_power.png'
    lv_active_power_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in lv_active_power_fsets:
        print(f"Saving {fset.title}")
        lv_active_power_recovery_plot(project_tuning, fset, data_root_dir, output_dirpath, plot_rise_time=True)




    # --- Healthy Phase Box Plots
    healthy_phase_fsets = []
    std_filter = ((analysis_df['Category'] != 'Fault_Ohms') & (analysis_df['Fault_Duration'] == 0.43))

    filter = (std_filter)
    title = "Healthy Phase Max. POC Voltage\n(All Faults)"
    filename = 'Healthy_Phase_Max_POC_V_pu_BoxPlot_AllFaults.png'
    healthy_phase_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    
    for fset in healthy_phase_fsets:
        print(f"Saving {fset.title}")
        healthy_phase_box_and_whisker(project_tuning, fset, data_root_dir, output_dirpath)

    # --- Max Voltage Post Fault Box Plots
    post_fault_V_phase_fsets = []
    std_filter = ((analysis_df['Category'] != 'Fault_Ohms') & (analysis_df['Fault_Duration'] == 0.43))
    
    filter = (std_filter)
    title = "Post-Fault Max. POC Voltage\n(All Faults)"
    filename = 'Post_Fault_Max_POC_V_pu_BoxPlot_AllFaults.png'
    post_fault_V_phase_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in post_fault_V_phase_fsets:
            print(f"Saving {fset.title}")
            max_phase_voltages_post_fault(project_tuning, fset, data_root_dir, output_dirpath)

    # --- Healthy Phase Scatter Plots
    healthy_phase_fsets = []
    filter = ((analysis_df['Init_Fault_MVA'] == max_flt_level) & (analysis_df['Fault_Duration'] == 0.43))
    title = "Healthy Phase Voltages\n(Maximum Fault Level [3138MVA])"
    filename = 'HealthyPhaseVoltages_ScatterPlot_MaxFaultLevel.png'
    healthy_phase_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = ((analysis_df['Init_Fault_MVA'] == min_flt_level) & (analysis_df['Fault_Duration'] == 0.43))
    title = "Healthy Phase Voltages\n(Minimum Fault Level [523MVA])"
    filename = 'HealthyPhaseVoltages_ScatterPlot_MinFaultLevel.png'
    healthy_phase_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in healthy_phase_fsets:
        print(f"Saving {fset.title}")
        healthy_phase_vs_fault_depth(project_tuning, fset, data_root_dir, output_dirpath)

    # --- Post Fault V Box plots
    post_fault_V_phase_fsets = []
    filter = (std_filter)
    title = "Post-Fault Max. POC Voltage\n(All Faults)"
    filename = 'Post_Fault_Max_POC_V_pu_BoxPlot_AllFaults.png'
    post_fault_V_phase_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in post_fault_V_phase_fsets:
        print(f"Saving {fset.title}")
        max_phase_voltages_post_fault(project_tuning, fset, data_root_dir, output_dirpath)
        
    # --- Settling Time Plots
    settling_time_fsets = []
    
    # - Non-resistive faults
    std_filter = ((analysis_df['Category'] != 'Fault_Ohms')  &
                  (analysis_df['Fault_Duration'] == 0.43) &
                  (analysis_df[FRT_MUST_ENGAGE_BY] == True) &
                  (pd.notna(analysis_df["Fault_iq_AEMO_Rise_Time"]))
                )       
    for fltlevel, fltlabel in [(min_flt_level, 'MinFaultLevel'), (max_flt_level, 'MaxFaultLevel')]:
        filter = (std_filter & (analysis_df['Init_Fault_MVA'] == fltlevel))
        title = f'Capacitive Current Injection\n(Non-Resistive low-voltage {fltlabel})'
        filename = f'Capacitive_iq_Non-Resistive{fltlabel}.png'
        settling_time_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
        
        filter = (std_filter &
                  (analysis_df['Init_Fault_MVA'] == fltlevel)
                  )
        title = f'Capacitive Current Injection\n(Slow Settling low-voltage {fltlabel})'
        filename = f'Capacitive_iq_SlowSettling_{fltlabel}.png'
        settling_time_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    
    # - Resistive faults    
    std_filter = ((analysis_df['Category'] == 'Fault_Ohms') &
                  (analysis_df['Fault_Duration'] == 0.43) &
                  (analysis_df[FRT_MUST_ENGAGE_BY] == True) &
                  (pd.notna(analysis_df["Fault_iq_AEMO_Rise_Time"]))
                  )
    for fltlevel, fltlabel in [(min_flt_level, 'MinFaultLevel'), (max_flt_level, 'MaxFaultLevel')]:
        filter = (std_filter & (analysis_df['Init_Fault_MVA'] == fltlevel))
        title = f'Capacitive Current Injection\n(Resistive low-voltage {fltlabel})'
        filename = f'Capacitive_iq_Resistive_{fltlabel}.png'
        settling_time_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
       
    # - All faults    
    std_filter = ((analysis_df['Fault_Duration'] == 0.43) &
                  (analysis_df[FRT_MUST_ENGAGE_BY] == True) &
                  (pd.notna(analysis_df["Fault_iq_AEMO_Rise_Time"]))
                  )
    for fltlevel, fltlabel in [(min_flt_level, 'MinFaultLevel'), (max_flt_level, 'MaxFaultLevel')]:
        filter = (std_filter & (analysis_df['Init_Fault_MVA'] == fltlevel))
        title = f'Capacitive Current Injection\n(All low-voltage {fltlabel})'
        filename = f'Capacitive_iq_All_{fltlabel}.png'
        settling_time_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in settling_time_fsets:
        print(f"Saving {fset.title}")
        settling_time_plot(project_tuning, fset, data_root_dir, output_dirpath, plot_rise_time=True)                      

    # --- Active Power Recovery Plots.
    active_power_recovery_fsets = []
    init_ppoc_mw = analysis_df['Init_Pwind_MW'] + analysis_df['Init_Pbess_MW']
    filter = ((analysis_df['Fault_Duration'] == 0.43) & (init_ppoc_mw == 400)
              & (analysis_df[FRT_MUST_ENGAGE_BY] == True) & (analysis_df['Xf_Ohms'].isnull()))
    title = f"Post-fault Active-Power Recovery\n Non-Resistive Faults Exporting (Init. Ppoc=400 MW)"
    filename = f'ActivePowerRecovery_Ures_Exporting.png'
    active_power_recovery_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = ((analysis_df['Fault_Duration'] == 0.43) & (init_ppoc_mw == -104) &
              (analysis_df[FRT_MUST_ENGAGE_BY] == True) & (analysis_df['Xf_Ohms'].isnull()))
    title = f"Post-fault Active-Power Recovery\n Non-Resistive Faults Importing (Init. Ppoc=-104 MW)"
    filename = f'ActivePowerRecovery_Ures_Importing.png'
    active_power_recovery_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = ((analysis_df['Fault_Duration'] == 0.43) & (init_ppoc_mw == 400)
              & (analysis_df[FRT_MUST_ENGAGE_BY] == True) & (analysis_df['Xf_Ohms'] == 0))
    title = f"Post-fault Active-Power Recovery\n Resistive Faults Exporting (Init. Ppoc=400 MW)"
    filename = f'ActivePowerRecovery_Resistive_Exporting.png'
    active_power_recovery_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = ((analysis_df['Fault_Duration'] == 0.43) & (init_ppoc_mw == -104) &
              (analysis_df[FRT_MUST_ENGAGE_BY] == True) & (analysis_df['Xf_Ohms'] == 0))
    title = f"Post-fault Active-Power Recovery\n Resistive Faults Importing (Init. Ppoc=-104 MW)"
    filename = f'ActivePowerRecovery_Resistive_Importing.png'
    active_power_recovery_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for gov_label, gov_thresh in [('no_freeze', analysis_df['Acc. GovFrz during fault'] == 0),
                                  ('little_freeze', analysis_df['Acc. GovFrz during fault'] < 0.2),
                                  ('fully_freeze', analysis_df['Acc. GovFrz during fault'] > 0.4)]:
        filter = ((analysis_df['Fault_Duration'] == 0.43) & (init_ppoc_mw == 400)
                  & (analysis_df[FRT_MUST_ENGAGE_BY] == True) & (analysis_df['Xf_Ohms'].isnull()))
        title = f"Post-fault Active-Power Recovery\n Non-Resistive Faults Exporting (Init. Ppoc=400 MW) [{gov_label}]"
        filename = f'ActivePowerRecovery_Ures_Exporting_{gov_label}.png'
        active_power_recovery_fsets.append(FigureSet(analysis_df[filter & gov_thresh], title, filename_prefix + filename))

        filter = ((analysis_df['Fault_Duration'] == 0.43) & (init_ppoc_mw == -104) &
                  (analysis_df[FRT_MUST_ENGAGE_BY] == True) & (analysis_df['Xf_Ohms'].isnull()))
        title = f"Post-fault Active-Power Recovery\n Non-Resistive Faults Importing (Init. Ppoc=-104 MW)[{gov_label}]"
        filename = f'ActivePowerRecovery_Ures_Importing_{gov_label}.png'
        active_power_recovery_fsets.append(FigureSet(analysis_df[filter & gov_thresh], title, filename_prefix + filename))

        filter = ((analysis_df['Fault_Duration'] == 0.43) & (init_ppoc_mw == 400)
                  & (analysis_df[FRT_MUST_ENGAGE_BY] == True) & (analysis_df['Xf_Ohms'] == 0))
        title = f"Post-fault Active-Power Recovery\n Resistive Faults Exporting (Init. Ppoc=400 MW) [{gov_label}]"
        filename = f'ActivePowerRecovery_Resistive_Exporting_{gov_label}.png'
        active_power_recovery_fsets.append(FigureSet(analysis_df[filter & gov_thresh], title, filename_prefix + filename))

        filter = ((analysis_df['Fault_Duration'] == 0.43) & (init_ppoc_mw == -104) &
                  (analysis_df[FRT_MUST_ENGAGE_BY] == True) & (analysis_df['Xf_Ohms'] == 0))
        title = f"Post-fault Active-Power Recovery\n Resistive Faults Importing (Init. Ppoc=-104 MW) [{gov_label}]"
        filename = f'ActivePowerRecovery_Resistive_Importing_{gov_label}.png'
        active_power_recovery_fsets.append(FigureSet(analysis_df[filter & gov_thresh], title, filename_prefix + filename))

    for fset in active_power_recovery_fsets:
        print(f"Saving {fset.title}")
        # Plot the Normal PLot.
        ppoc_error_settling_time_plot(project_tuning, fset, data_root_dir, output_dirpath)

    # --- Settling Time Plots
    rise_and_settle_against_vpoc_deviation = []

    for FRT_ENGAGE in ['Fault_Vpoc_below_80', 'Fault_Vpoc_below_85']:
        filter = ((analysis_df['Fault_Duration'] == 0.43) & (analysis_df[FRT_ENGAGE] == True))
        title = f'All Faults'
        filename = f'Rise_and_Settle_Voltage_Deviation_All_Faults_{FRT_ENGAGE}.png'
        rise_and_settle_against_vpoc_deviation.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for FRT_ENGAGE in ['Fault_Vpoc_below_80', 'Fault_Vpoc_below_85']:
        filter = ((analysis_df['Fault_Duration'] == 0.43) & (analysis_df[FRT_ENGAGE] == True)
                  & ~(analysis_df['Category'] == 'Fault_Ohms'))
        title = f'Non-Resitive Faults'
        filename = f'Rise_and_Settle_Voltage_Deviation_All_Faults_{FRT_ENGAGE}.png'
        rise_and_settle_against_vpoc_deviation.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in rise_and_settle_against_vpoc_deviation:
        print(f"Saving {fset.title}")
        # Plot the Normal PLot.
        rise_and_settle_times_against_fault_delta_v(project_tuning, fset, data_root_dir, output_dirpath)

if __name__ == '__main__':
    standard_plotting_main(plotting_function=plotting_function)

