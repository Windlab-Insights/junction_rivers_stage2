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
        raise ValueError("Inductive Iq Injection Settling time plot requires figure set to pass in sets of the same Fault_Duration.")
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
                ]
        blabels = ['Settling-time', 'Rise-time']
    else:
        bdata = [list(df_subset[f'Fault_iq_AEMO_Settling_Time'].values),
                 ]
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
    fig.set_size_inches(12 * cm, 12 * cm)
    fig.subplots_adjust(hspace=0, left=0.18, right=0.95, bottom=0.1)

    max_sig_settling_time = max(df_subset[f'Ppoc_95pc_Recovery_Time'])

    ax = axs[0]
    axs[0].grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)

    color = iter(plt.cm.inferno(np.linspace(0, 0.8, len(df_subset))))
    for _, scenario in df_subset.sort_values(by=['Fault_Vpoc_pu'], ascending=False).iterrows():
        step_start = model_init_time + scenario['Fault_Time']
        step_end = step_start + fault_duration

        plot_start = step_end - (1 / 8) * max_sig_settling_time
        plot_end = plot_start + max((2.2) * max_sig_settling_time, 0.1)
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

    # x_vals = np.linspace(plot_start- step_end, plot_end - step_end, 100)
    # settle_line = -0.05 * ((scenario['Init_Pwind_MW'] + scenario['Init_Pbess_MW']) / ppoc_mw_base)
    # y_vals = [settle_line] * len(x_vals)
    # fill_vals = [np.sign(settle_line)] * len(x_vals)
    # axs[0].axhline(settle_line, color='black', linestyle='--', linewidth='0.5')
    # axs[0].axhline(0, color='black', linestyle='--', linewidth='0.5')
    # axs[0].fill_between(x_vals, y_vals, fill_vals, color='yellow', alpha=0.05)

    axs[0].set_ylim(curr_ylims)
    new_metric = np.minimum(df_subset[f'Ppoc_Within10MW_Recovery_Time'].values, df_subset[f'Ppoc_95pc_Recovery_Time'].values)
    bdata = [list(df_subset[f'Ppoc_95pc_Recovery_Time'].values)]
    blabels = ['95%\nrecovery']
    # bdata = [list(new_metric), list(df_subset[f'Ppoc_95pc_Recovery_Time'].values)]
    # blabels = ['95% or 10MW', '95%\nrecovery']

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


def high_voltage_iq_injection_plot_plq_profile(tuning, figure_set, pkl_dir, output_dir):

    plq_profile_spec_path = r"D:\grid_workspace\gawarabaya\specs\specs__csr_s5254_plq_voltage_profile.csv"
    plq_profile_spec_df = pd.read_csv(plq_profile_spec_path)
    plq_profile_pkl_dir = r"G:\GBWF_PSCAD_Results\2024-01-19-PSCAD-GBWF-v3.25-results\CSR_S5254"
    
    plq_profile_vpoc_values = []
    plq_profile_delta_iq_values = []

    for _, scenario in plq_profile_spec_df.iterrows():
        
        plq_profile_pkl_path = get_file_path(plq_profile_pkl_dir,  scenario["File_Name"] + ".pkl")

        plq_profile_df = pd.read_pickle(plq_profile_pkl_path)

        for poc_v, poc_iq in zip(plq_profile_df["POC_V_pu"].values,  plq_profile_df["POC_Iq_pos_pu"].values):

            print(poc_v, poc_iq)

            if poc_v == 1.25 or poc_v == 1.30:

                print("-----",poc_v, poc_iq)

                plq_profile_vpoc_values.append(poc_v)
                plq_profile_delta_iq_values.append(poc_iq)  

                break




    df_subset = figure_set.scenarios
    if len(df_subset) == 0:
        return

    init_matplotlib()
    fig_title = figure_set.title
    plot_filename = figure_set.filename

    def hvrt_iq_guide(k, u, v):
        return np.clip(k*(u - v), -1, 0)

    v_array = np.linspace(1, 1.5, 51*6)
    iq_mas_values = [hvrt_iq_guide(0, 1.2, v) for v in v_array]
    iq_mas_old_values = [hvrt_iq_guide(3, 1.15, v) for v in v_array]
    iq_ass_values = [hvrt_iq_guide(6, 1.15, v) for v in v_array]

    figs, ax = plt.subplots(1)
    cm = 1 / 2.54
    figs.set_size_inches(12 * cm, 12 * cm)
    figs.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)

    ax.plot(v_array, iq_mas_values, ls='-', color='yellow', label='mas[0%/%]')
    ax.plot(v_array, iq_mas_old_values, ls='-', color='gray', label='nas[3%/%]')
    ax.plot(v_array, iq_ass_values, ls='-', color='green', label='aas[6%/%]')

    ax.fill_between(v_array, iq_mas_old_values, iq_mas_values, color='yellow', alpha=0.1)
    ax.fill_between(v_array, iq_ass_values, iq_mas_old_values, color='gray', alpha=0.1)
    ax.fill_between(v_array, [-2] * len(v_array), iq_ass_values, color='green', alpha=0.1)
    ax.grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)

    ax.scatter(df_subset['Fault_Vpoc_pu'], df_subset['Delta_iq_mcc_pu'], edgecolor='black', facecolor='white', zorder=2)
    ax.scatter(plq_profile_vpoc_values, plq_profile_delta_iq_values, edgecolor='blue', facecolor='white', zorder=2)

    ax.set_title(fig_title)
    ax.set_ylabel(f"$i_q$ Injection [pu mcc]")
    ax.set_xlabel("POC Settled Fault Voltage [pu]")
    ax.set_xlim(1.0, 1.4)
    ax.set_ylim(-1.2, 0.1)
    ax.legend()
    figs.savefig(os.path.join(output_dir, plot_filename))




def high_voltage_iq_injection_plot(tuning, figure_set, pkl_dir, output_dir):
    df_subset = figure_set.scenarios
    if len(df_subset) == 0:
        return

    init_matplotlib()
    fig_title = figure_set.title
    plot_filename = figure_set.filename

    def hvrt_iq_guide(k, u, v):
        return np.clip(k*(u - v), -1, 0)

    v_array = np.linspace(1, 1.5, 51*6)
    iq_mas_values = [hvrt_iq_guide(0, 1.2, v) for v in v_array]
    iq_mas_old_values = [hvrt_iq_guide(3, 1.15, v) for v in v_array]
    iq_ass_values = [hvrt_iq_guide(6, 1.15, v) for v in v_array]

    figs, ax = plt.subplots(1)
    cm = 1 / 2.54
    figs.set_size_inches(12 * cm, 12 * cm)
    figs.subplots_adjust(hspace=0, left=0.15, right=0.95, bottom=0.1)

    ax.plot(v_array, iq_mas_values, ls='-', color='yellow', label='mas[0%/%]')
    ax.plot(v_array, iq_mas_old_values, ls='-', color='gray', label='nas[3%/%]')
    ax.plot(v_array, iq_ass_values, ls='-', color='green', label='aas[6%/%]')

    ax.fill_between(v_array, iq_mas_old_values, iq_mas_values, color='yellow', alpha=0.1)
    ax.fill_between(v_array, iq_ass_values, iq_mas_old_values, color='gray', alpha=0.1)
    ax.fill_between(v_array, [-2] * len(v_array), iq_ass_values, color='green', alpha=0.1)
    ax.grid(which='major', axis='both', linestyle='--', color='lightgray', alpha=0.6)

    ax.scatter(df_subset['Fault_Vpoc_pu'], df_subset['Delta_iq_mcc_pu'], edgecolor='black', facecolor='white', zorder=2)
    ax.set_title(fig_title)
    ax.set_ylabel(f"$i_q$ Injection [pu mcc]")
    ax.set_xlabel("POC Settled Fault Voltage [pu]")
    ax.set_xlim(1.0, 1.4)
    ax.set_ylim(-1.2, 0.1)
    ax.legend()
    figs.savefig(os.path.join(output_dir, plot_filename))



def plotting_function(project_tuning, analysis_df, data_root_dir, output_dirpath):

    filename_prefix = "s5255_high_voltage_"

    min_flt_level = 523
    max_flt_level = 3138
    FRT_MUST_ENGAGE_BY = 'Fault_Vpoc_above_115'

    # --- high Voltage IQ Injection plots
    iq_injection_fsets = []
    std_filter = ((analysis_df['Fault_Duration'] == 0.43) &
                  (analysis_df['Fault_iq_mcc_pu'] > -1.0) &
                  (analysis_df[FRT_MUST_ENGAGE_BY] == True))
    filter = (std_filter & (analysis_df['Init_Fault_MVA'] == max_flt_level))
    title = "POC $i_q$ Injection: Max Fault Level [3138MVA]\nOver-voltage Faults"
    filename = 'POC_Iq_Injection_MaxFaultLevel.png'
    iq_injection_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = (std_filter & (analysis_df['Init_Fault_MVA'] == min_flt_level))
    title = "POC $i_q$ Injection: Min. Fault Level [523MVA]\nOver-voltage Faults"
    filename = 'POC_Iq_Injection_MinFaultLevel.png'
    iq_injection_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))


    title = "POC $i_q$ Injection: Over-voltage Faults"
    filename = 'POC_Iq_Injection_all.png'
    iq_injection_fsets.append(FigureSet(analysis_df, title, filename_prefix + filename))

    high_voltage_iq_injection_plot_plq_profile(project_tuning, FigureSet(analysis_df, title, filename_prefix + filename), data_root_dir, output_dirpath)

    for fset in iq_injection_fsets:
        print(f"Saving {fset.title}")
        high_voltage_iq_injection_plot(project_tuning, fset, data_root_dir, output_dirpath)

    



    
    # --- Settling Time Plots
    settling_time_fsets = []
    for fltlevel,fltlabel in [(min_flt_level, 'MinFaultLevel'), (max_flt_level, 'MaxFaultLevel')]:
        filter = ((analysis_df['Fault_Duration'] == 0.43) &
                (analysis_df[FRT_MUST_ENGAGE_BY] == True) & (analysis_df['Init_Fault_MVA'] == fltlevel) &
                (analysis_df[f'Fault_iq_AEMO_Settling_Time'].notna()))
        title = f'Inductive Current Injection\n(All faults high-voltage {fltlabel})'
        filename = f'Inductive_iq_All_{fltlabel}.png'
        settling_time_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    
    filter = ((analysis_df['Fault_Duration'] == 0.43) &
            (analysis_df[FRT_MUST_ENGAGE_BY] == True) & 
            (analysis_df[f'Fault_iq_AEMO_Settling_Time'].notna()))
    title = f'Inductive Current Injection\n(All faults high-voltage)'
    filename = f'Inductive_iq_All.png'
    settling_time_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    
    
    for fset in settling_time_fsets:
        print(f"Saving {fset.title}")
        # Plot the Normal PLot.
        settling_time_plot(project_tuning, fset, data_root_dir, output_dirpath, plot_rise_time=True)

    # --- Active Power Recovery Plots.
    active_power_recovery_fsets = []

    init_ppoc_mw = analysis_df['Init_Pwind_MW'] + analysis_df['Init_Pbess_MW']
    filter = ((analysis_df['Fault_Duration'] == 0.43)
              & (analysis_df[FRT_MUST_ENGAGE_BY] == True))
    title = f"Post Over-Voltage Active-Power Recovery"
    filename = f'ActivePowerRecovery.png'
    active_power_recovery_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    init_ppoc_mw = analysis_df['Init_Pwind_MW'] + analysis_df['Init_Pbess_MW']
    filter = ((analysis_df['Fault_Duration'] == 0.43) & (init_ppoc_mw == 400)
              & (analysis_df[FRT_MUST_ENGAGE_BY] == True))
    # print(filter)
    title = f"Post Over-Voltage Active-Power Recovery\n (Init. Ppoc=400 MW)"
    filename = f'ActivePowerRecovery_Exporting.png'
    active_power_recovery_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))


    init_ppoc_mw = analysis_df['Init_Pwind_MW'] + analysis_df['Init_Pbess_MW']
    filter = ((analysis_df['Fault_Duration'] == 0.43) & (init_ppoc_mw == -104) &
              (analysis_df[FRT_MUST_ENGAGE_BY] == True))
    title = f"Post Over-Voltage Active-Power Recovery\n (Init. Ppoc=-104 MW)"
    filename = f'ActivePowerRecovery_Importing.png'
    active_power_recovery_fsets.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in active_power_recovery_fsets:
        print(f"Saving {fset.title}")
        ppoc_error_settling_time_plot(project_tuning, fset, data_root_dir, output_dirpath)


if __name__ == '__main__':
    standard_plotting_main(plotting_function=plotting_function)

