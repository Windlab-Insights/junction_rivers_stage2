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

def plot_qpoc_vpoc_grid(tuning, figure_set: FigureSet, pkl_dir, output_dir):
    df_subset = figure_set.scenarios
    fig_title = figure_set.title
    plot_filename = figure_set.filename

    # print("figure_set.scenarios: ",figure_set.scenarios)

    init_matplotlib()

    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])

    nrows = 5
    ncols = 3
    total_plots = len(df_subset)
    plots_per_page = (nrows * ncols)
    total_plot_pages = math.ceil(len(df_subset) / plots_per_page)
    plot_pages = []
    for p_idx in range(total_plot_pages):
        fig, axs = plt.subplots(nrows, ncols)
        plt.subplots_adjust(left=0.05,
                            bottom=0.08,
                            right=0.95,
                            top=0.92,
                            wspace=0.45,
                            hspace=0.43)
        fig.set_size_inches(11.7, 8.3)
        plot_pages.append((fig, axs))

    plot_id = -1
    for _, scenario in df_subset.iterrows():
        plot_id += 1
        curr_plot_page = plot_id // plots_per_page
        curr_plot_id = plot_id % plots_per_page

        fig, axs = plot_pages[curr_plot_page]
        ax = axs.reshape(-1)[curr_plot_id]

        osc_freq = float(scenario.Osc_Freq_Hz)
        fault_time = float(scenario.Fault_Time)
        fault_duration = float(scenario.Fault_Duration)

        pkl_path = get_file_path(pkl_dir, scenario.File_Name + ".pkl")

        if os.path.isfile(pkl_path):
            data = pd.read_pickle(pkl_path)
            osc_start_time = model_init_time + fault_time
            osc_end_time = osc_start_time + fault_duration
            plot_length = min((3 / osc_freq), osc_end_time - osc_start_time)
            plot_start_time = osc_start_time - plot_length / 8
            plot_end_time = osc_start_time + plot_length

            v_signal = data["POC_V_pu"]
            q_signal = data["POC_Q_MVAr"]

            ax.plot(v_signal[plot_start_time:plot_end_time].index, v_signal[plot_start_time:plot_end_time].values,
                    color='r')
            ax2 = ax.twinx()
            ax2.plot(q_signal[plot_start_time:plot_end_time].index, q_signal[plot_start_time:plot_end_time].values,
                     '--', color='b')

            minx = plot_start_time
            # maxy = max(v_signal[plot_start_time:plot_end_time].values)
            maxy = ax.get_ylim()[1]
            ax.text(minx, maxy, f'{osc_freq:.1f} hz', ha='left', va='top', backgroundcolor='white')
            ax.set_xlim(plot_start_time, plot_end_time)
        else:
            continue
            raise Exception(f"No File found={pkl_path}")

    for page_idx, (fig, axs) in enumerate(plot_pages):
        fig.suptitle(f'Oscillation disturbance Response of Vpoc (left-axis red-solid .pu) and Qpoc (right-axis blue-dashed MVAr).\n{fig_title}')
        fig_name = f"{plot_filename}__{page_idx+1}of{len(plot_pages)}.png"
        outfile = os.path.join(output_dir, fig_name)
        fig.savefig(outfile)
        plt.close(fig)

# ------------


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]


def generate_freq_response(tuning, df_subset, pkl_dir):
    vslack_to_vpoc_freqs = []
    vslack_to_vpoc_mag = []
    vpoc_to_qpoc_phase = []
    for _, scenario in df_subset.iterrows():
        # scenario = df_subset.loc[i]
        osc_freq = scenario.Osc_Freq_Hz
        osc_amp = scenario.Osc_Amplitude
        fault_time = scenario.Fault_Time
        fault_duration = scenario.Fault_Duration
        model_init_time = float(tuning['TIME_Full_Init_Time_sec'])

        pkl_path = get_file_path(pkl_dir, scenario.File_Name + ".pkl")

        print(f"Pickle File: {pkl_path}")
    

        if pkl_path == False:
            continue

        data = pd.read_pickle(pkl_path)

        osc_start_time = model_init_time + fault_time
        osc_end_time = osc_start_time + fault_duration

        v_signal = data["POC_V_pu"][osc_start_time:osc_end_time]
        q_signal = data["POC_Q_MVAr"][osc_start_time:osc_end_time]

        sample_period = v_signal.index[1] - v_signal.index[0]
        v_fft_spectrum = scipy.fft.rfft(v_signal.values)
        v_fft_freq = scipy.fft.rfftfreq(len(v_signal), sample_period)
        q_fft_spectrum = scipy.fft.rfft(q_signal.values)
        q_fft_freq = scipy.fft.rfftfreq(len(q_signal), sample_period)

        idx, val = find_nearest(v_fft_freq, osc_freq)
        vslack_to_vpoc_freqs.append(val)
        v_norm_coeff = 2 / (len(v_signal) * osc_amp)
        vslack_to_vpoc_mag.append(abs(v_fft_spectrum[idx] * v_norm_coeff))

        phase_diff = np.rad2deg(cmath.phase(q_fft_spectrum[idx]) - cmath.phase(v_fft_spectrum[idx]))
        vpoc_to_qpoc_phase.append(phase_diff)

    return vslack_to_vpoc_freqs, vslack_to_vpoc_mag, vslack_to_vpoc_freqs, vpoc_to_qpoc_phase




def plot_freq_response(tuning, figure_set: FigureSet, pkl_dir, output_dir):
    df_subset = figure_set.scenarios
    fig_title = figure_set.title
    plot_filename = figure_set.filename

    # print("\nPlot Freq Response: \n\tfigure_set.scenarios: ",figure_set.scenarios,"\n")

    init_matplotlib()
    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])

    fig, axs = plt.subplots(2, 1)
    output = generate_freq_response(tuning, df_subset, pkl_dir)
    vslack_to_vpoc_freqs, vslack_to_vpoc_mag, _, vpoc_to_qpoc_phase = output

    axs[0].axhline(0, linestyle='--', color='black', alpha=0.8, linewidth=0.7)
    # axs[0].add_patch(matplotlib.patches.Rectangle((0, -100), 1000, 100, alpha=0.1, color='lightgreen'))
    # p = matplotlib.patches.Rectangle((0, 0), 1000, 100, alpha=0.3, color='red')
    # axs[0].add_patch(p)

    axs[0].scatter(vslack_to_vpoc_freqs, 20 * np.log10(vslack_to_vpoc_mag), marker='o', s=3, color='black')
    axs[0].set_title(fig_title)
    axs[0].set_xscale('log')
    axs[0].grid(which='both', linestyle='--', linewidth=0.7)
    axs[0].set(xlim=(0.07, 50), ylim=(-50, 10))
    axs[0].set_ylabel(r'Gain Vgrid $\rightarrow$ Vpoc (dB)')

    axs[1].grid(which='both', linestyle='--', linewidth=0.7, zorder=0)
    # axs[1].axhline(-180, linestyle=':', color='g', alpha=0.8, linewidth=1, zorder=1)
    # axs[1].scatter(vslack_to_vpoc_freqs, np.unwrap(vpoc_to_qpoc_phase, 360), marker=3)
    axs[1].scatter(vslack_to_vpoc_freqs, vpoc_to_qpoc_phase, marker='o', s=3, color='black')
    ones_array = np.array([1] * len(vslack_to_vpoc_freqs))
    axs[1].fill_between(axs[0].get_xlim(), [90, 90], [180, 180], alpha=0.2, color='green')
    axs[1].fill_between(axs[0].get_xlim(), [-90, -90], [-180, -180], alpha=0.2, color='green')
    axs[1].set_xscale('log')
    axs[1].set(xlim=(0.07, 50), ylim=(-180, 180))
    axs[1].set_ylabel(r'ph(Qpoc) - ph(Vpoc) (deg)')
    axs[1].set_yticks([-180, -90, 0, 90, 180])
    plt.xlabel(r'Frequency ($\rm hz$)')

    fig.align_labels()
    outfile = os.path.join(output_dir, plot_filename)
    plt.savefig(outfile)
    plt.close()



def plotting_function(project_tuning, analysis_df, data_root_dir, output_dirpath):

    # print("\nanalysis_df: ",analysis_df)

    filename_prefix = "osc rej "

    min_flt_level = 523
    max_flt_level = 3138

    qpoc_vpoc_grid_plots = []
    filter = ((analysis_df['Osc_Phase_deg'] == 0) & (analysis_df['Init_Fault_MVA'] == max_flt_level))
    title = '[Max. Fault Level]'
    filename = 'VpocQpocGridPlot_MaxFaultLvl'
    qpoc_vpoc_grid_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = ((analysis_df['Osc_Phase_deg'] == 0) & (analysis_df['Init_Fault_MVA'] == min_flt_level))
    title = '[Min. Fault Level]'
    filename = 'VpocQpocGridPlot_MinFaultLvl'
    qpoc_vpoc_grid_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in qpoc_vpoc_grid_plots:
        print(f"Saving {fset.title}")
        plot_qpoc_vpoc_grid(project_tuning, fset, data_root_dir, output_dirpath)

    freq_resp_plots = []
    filter = ((analysis_df['Init_Fault_MVA'] == max_flt_level) & (analysis_df['Osc_Phase_deg'] == 0))
    title = 'Vgrid to POC Gain & Phase Analysis\nMax. Fault Level '
    filename = 'OscRej_FreqResp_MaxFaultLevel.png'
    freq_resp_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    
    filter = ((analysis_df['Init_Fault_MVA'] == min_flt_level) & (analysis_df['Osc_Phase_deg'] == 0))
    title = 'Vgrid to POC Gain & Phase Analysis\Min. Fault Level '
    filename = 'OscRej_FreqResp_MaxFaultLevel_2deg.png'
    freq_resp_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    filter = ((analysis_df['Init_Fault_MVA'] == max_flt_level) & (analysis_df['Osc_Phase_deg'] == 2))
    title = 'Vgrid to POC Gain & Phase Analysis\nMax. Fault Level / 2deg '
    filename = 'OscRej_FreqResp_MinFaultLevel.png'
    freq_resp_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    
    filter = ((analysis_df['Init_Fault_MVA'] == min_flt_level) & (analysis_df['Osc_Phase_deg'] == 2))
    title = 'Vgrid to POC Gain & Phase Analysis\Min. Fault Level / 2deg '
    filename = 'OscRej_FreqResp_MinFaultLevel_2deg.png'
    freq_resp_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    #
    # filter = ((analysis_df['Init_Fault_MVA'] == min_flt_level) & (analysis_df['Init_Pwind_MW'] == 400) & (analysis_df['Osc_Phase_deg'] == 0))
    # title = 'Vgrid to POC Gain & Phase Analysis\nMax. Fault Level | Pwind(400) | Pbess(0)'
    # filename = 'OscRej_FreqResp_MinFaultLevel_Pwind_400_Pbess_0.png'
    # freq_resp_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    
    # filter = ((analysis_df['Init_Fault_MVA'] == max_flt_level) & (analysis_df['Init_Pwind_MW'] == 400) & (analysis_df['Osc_Phase_deg'] == 2))
    # title = 'Vgrid to POC Gain & Phase Analysis\nMax. Fault Level | Pwind(400) | Pbess(0) | 2deg'
    # filename = 'OscRej_FreqResp_MaxFaultLevel_Pwind_400_Pbess_0_2deg.png'
    # freq_resp_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    
    # filter = ((analysis_df['Init_Fault_MVA'] == min_flt_level) & (analysis_df['Init_Pwind_MW'] == 400) & (analysis_df['Osc_Phase_deg'] == 2))
    # title = 'Vgrid to POC Gain & Phase Analysis\nMin. Fault Level | Pwind(400) | Pbess(0) | 2deg'
    # filename = 'OscRej_FreqResp_MinFaultLevel_Pwind_400_Pbess_0_2deg.png'
    # freq_resp_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))

    for fset in freq_resp_plots:
        print(f"Saving {fset.title}")
        plot_freq_response(project_tuning, fset, data_root_dir, output_dirpath)


    # ----- ----------
    #
    # # --- Ppoc Plots
    # ppoc_settling_plots = []
    # # exclude_categories = "|".join([re.escape(x) for x in ["Cont.Ctrl", "No-Windup[D.wc]", "No-Windup[U.wc]"]])
    # # filter = (~analysis_df['Category'].str.contains(exclude_categories))
    # filter = (~analysis_df['Category'].isin(["Cont.Ctrl", "No-Windup[D.wc]", "No-Windup[U.wc]"]))
    # title = 'Vref Step: Ppoc Error\n(All Tests)'
    # filename = 'VrefSteps_All_Ppoc_Error.png'
    # ppoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    #
    # # include_categories = "|".join([re.escape(x) for x in ["Sat(d5%)", "Sat(u5%)"]])
    # # filter = (analysis_df['Category'].str.contains(include_categories))
    # filter = (analysis_df['Category'].isin(["d5%", "u5%"]))
    # title = f'Vref Step: Ppoc Error\n(5%-step No Saturation)'
    # filename = 'VrefSteps_5pc_NoSaturation_Ppoc_Error.png'
    # ppoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    #
    # # include_categories = "|".join([re.escape(x) for x in ["PLQSat(d5%)", "PLQSat(u5%)"]])
    # # filter = (analysis_df['Category'].str.contains(include_categories))
    # filter = (analysis_df['Category'].isin(["Sat(d5%)", "Sat(u5%)"]))
    # title = f'Vref Step: Ppoc Error\n(5%-step Saturation)'
    # filename = 'VrefSteps_5pc_Saturation_Ppoc_Error.png'
    # ppoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    #
    #
    # for fset in ppoc_settling_plots:
    #     print(f"Saving {fset.title}")
    #     vref_ppoc_error_plot(project_tuning, fset, data_root_dir, output_dirpath)
    #
    #
    # # --- Vpoc Plots
    # vpoc_settling_plots = []
    # # exclude_categories = "|".join([re.escape(x) for x in ["Cont.Ctrl", "No-Windup[D.wc]", "No-Windup[U.wc]"]])
    # # filter = (~analysis_df['Category'].str.contains(exclude_categories))
    # filter = (~analysis_df['Category'].isin(["Cont.Ctrl", "No-Windup[D.wc]", "No-Windup[U.wc]"]))
    # title = 'Vref Step: Vpoc\n(All Tests)'
    # filename = 'VrefSteps_All_Vpoc.png'
    # vpoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    #
    # # filter = (analysis_df['Category'].str.contains(include_categories))
    # filter = (analysis_df['Category'].isin(["d5%", "u5%"]))
    # title = 'Vref Step: Vpoc\n(5%-step No Saturation)'
    # filename = 'VrefSteps_5pc_No_Saturation_Vpoc.png'
    # vpoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    #
    # include_categories = "|".join([re.escape(x) for x in ["Sat(d5%)", "Sat(u5%)"]])
    # filter = (analysis_df['Category'].isin(["Sat(d5%)", "Sat(u5%)"]))
    # title = f'Vref Step: Vpoc\n(5%-step Saturation)'
    # filename = 'VrefSteps_5pc_Saturation_Vpoc.png'
    # vpoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    #
    # for fset in vpoc_settling_plots:
    #     print(f"Saving {fset.title}")
    #     vref_vpoc_plot(project_tuning, fset, data_root_dir, output_dirpath)
    #
    #
    # # --- Qpoc Plots
    # qpoc_settling_plots = []
    # exclude_categories = ["Cont.Ctrl", "No-Windup[D.wc]", "No-Windup[U.wc]"]
    # filter = (~analysis_df['Category'].isin(exclude_categories))
    # title = 'Vref Step: Qpoc Error\n(All Tests)'
    # filename = 'VrefSteps_All_Qpoc_Error.png'
    # qpoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    #
    # # Naming of Sat here is confusing, it's based on mark/windlab definition. They now correspond to tests that
    # # don't saturate according to powerlink/TNSP definition.
    # # include_categories = "|".join([re.escape(x) for x in ["Sat(d5%)", "Sat(u5%)"]])
    # include_categories = ["d5%", "u5%"]
    # filter = (analysis_df['Category'].isin(include_categories))
    # title = 'Vref Step: Qpoc Error\n(5%-step No Saturation)'
    # filename = 'VrefSteps_5pc_No_Saturation_Qpoc_Error.png'
    # qpoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    #
    #
    # # include_categories = "|".join([re.escape(x) for x in ["PLQSat(d5%)", "PLQSat(u5%)"]])
    # include_categories = ["Sat(d5%)", "Sat(u5%)"]
    # filter = (analysis_df['Category'].isin(include_categories))
    # title = f'Vref Step: Qpoc Error\n(5% Step Saturation)'
    # filename = 'VrefSteps_5pc_Saturation_Qpoc_Error.png'
    # qpoc_settling_plots.append(FigureSet(analysis_df[filter], title, filename_prefix + filename))
    #
    # for fset in qpoc_settling_plots:
    #     print(f"Saving {fset.title}")
    #     vref_qpoc_error_plot(project_tuning, fset, data_root_dir, output_dirpath)

if __name__ == '__main__':
    standard_plotting_main(plotting_function=plotting_function)

