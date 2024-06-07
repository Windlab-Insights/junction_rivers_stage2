import math
from icecream import ic
from rengen.plotting.Plotter import Plotter
import os
import json
import logging
from typing import List, Optional, Union, Dict
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from matplotlib.table import Cell, Table
from scipy.signal import butter, lfilter
from enum import Enum, auto

plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
plt.rcParams['axes.xmargin'] = 0
plt.rcParams['axes.labelsize'] = 'medium'
plt.rcParams['axes.titlesize'] = 'medium'
plt.rcParams['axes.grid'] = True
plt.rcParams['axes.grid.axis'] = 'both'
plt.rcParams['figure.titleweight'] = 'bold'
plt.rcParams['figure.dpi'] = 300

DECIMATE = 50
DPI = 300

AXIS_COLOUR = (0, 0, 0)
PRIMARY_COLOUR = (1, 1, 1)
SECONDARY_COLOUR = (0.98, 0.98, 0.98)

COL_REF = (0.87, 0.70, 0.06, 0.8)
COL_POC = (0.25, 0.33, 0.83, 0.8)
COL_SIG_1 = (0.71, 0.11, 0.08, 0.8)
COL_SIG_2 = (0.00, 0.75, 1.00, 0.8)
COL_SIG_3 = (0.00, 0.70, 0.37, 0.8)
COL_SIG_4 = (0.98, 0.29, 0.69, 0.8)
COL_FRZ = (0.6, 0.6, 1.0, 0.8)
COL_HVRT = (1.0, 0.2, 0.2, 0.8)
COL_LVRT = (0.2, 0.2, 1.0, 0.8)
COL_TRIP = (0.0, 0.0, 0.0, 0.8) 
COL_CENTER_LINE = (0.8, 0.8, 0.8, 0.8)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("UbwfAPscadPlotterV2")
logger.setLevel(logging.INFO)


class BBWFPlotterV1(Plotter):

    def plot_from_df_and_dict(
            self,
            df: pd.DataFrame,
            spec_dict: Dict,
            png_path: Optional[Union[Path, str]] = None,
            pdf_path: Optional[Union[Path, str]] = None,
    ):
        plt.clf()
        
        SIGNAL_MAP = [
            "POC_P_MW", 
            "POC_WTG_Pref_MW",
            "POC_Q_MVAr", 
            "Qref_droop_MVAr",
            "POC_Vrms_pu", 
            "POC_Vref_pu",
            "poc_freq_hz",
        ]
        
        init_time = 0 #float(spec_dict["substitutions"]["TIME_Full_Init_Time_sec"])
        plot_duration = df[SIGNAL_MAP[0]].index[-1] - init_time

        

        num_plots = len(SIGNAL_MAP)
        number_of_cols = 4
        number_of_rows = 6
        number_of_pages = math.ceil(num_plots/(number_of_rows*number_of_cols))

        for page_number in range(1,number_of_pages+1):

            plt.cla()
            plt.clf()

            signal_map = SIGNAL_MAP[(page_number-1)*(number_of_rows*number_of_cols):(page_number)*(number_of_rows*number_of_cols)]

            print(f"page {page_number} of {number_of_pages}")

            fig = mpl_setup(60,40)

            gs_main = gridspec.GridSpec(number_of_rows+1, number_of_cols, figure=fig)

            axis_list: List[plt.Axes] = []

            for col_index in range(number_of_cols):
                for row_index in range(number_of_rows):
                    axis_list.append(fig.add_subplot(gs_main[row_index+1,col_index]))
                    if (row_index+number_of_rows*col_index) > (len(signal_map)-1):
                        axis_list[-1].axis("off")

            # Plot Title    
            ax_title = fig.add_subplot(gs_main[0])
            ax_title.axis("off")
            ax_title.set_title(f" ", style = 'italic', fontsize='x-small', loc='center', y=1.00, pad=5, color=AXIS_COLOUR)

            # Table Data
            table_data = [
                ['Page', f"{page_number} of {number_of_pages}"],
                ['Filename/Int. Ref.', f"{spec_dict['File_Name']}"],
                ['Fault Level [MVA]', f"{spec_dict['Init_Fault_MVA']}"],
                ['X2R', f"{spec_dict['Init_Fault_X_on_R']}"],
                ['Vpoc [.pu]', f"{spec_dict['Init_Vpoc_pu']:.4f}"],
                ['Qpoc [.pu]', f"{spec_dict['Init_Qpoc_pu']}"],
                ['Pwind [MW]', f"{spec_dict['Init_Pwind_MW']}"],
                ['Pbess [MW]', f"{spec_dict['Init_Pbess_MW']}"],
            ]

            # Optional Table Data
            optional_table_information = [
                ('Vslack_pu', 'Vslack [.pu]'),
                ('Category', 'Category'),
                ('Fault_Type', 'Fault Type'),
                ('Fault_Time', 'Fault Time'),
                ('Fault_Duration', 'Fault Duration'),
                ('Fault_Zs_Multiplier', 'Zf/Zs'),
                ('Rf_Ohms', 'Rf'),
                ('Xf_Ohms', 'Xf'),
                ('Test_Profile', 'Test Profile'),
            ]

            for col_name, label in optional_table_information:
                if col_name in spec_dict:
                    table_data.append([label, f"{spec_dict[col_name]}"])

            # Plot Table
            table_colours = np.empty_like(table_data, dtype='object')
            for i,_ in enumerate(table_colours):
                table_colours[i] = [SECONDARY_COLOUR,SECONDARY_COLOUR]

            colWidths = [0.30, 0.40]
            ax_title.table(cellText=table_data, fontsize='x-small', colWidths=colWidths,cellLoc='left', loc='upper left',cellColours=table_colours)


            for signal_name, axis in zip(signal_map,axis_list):
                axis.cla()
                # axis = fig.add_subplot(gs_main[i])

                
                generate_std_subplot(axis, [signal_name], df, plot_duration, init_time, plot_start=0)
                axis.set_title(f"{signal_name}", style = 'italic', fontsize='x-small', loc='center' ,y=1.00, pad=5, color=AXIS_COLOUR)
            
            save_plot(png_path + f"_page_{page_number}", pdf=True)
            


def wrap_angle(signal):

    new_signal = signal

    for i, value in enumerate(signal):

        angle = value % 360
        angle = (angle + 360) % 360

        if (angle > 180):
            angle -= 360

        new_signal[i] = angle

    return signal



def generate_std_subplot(ax, signal_map, dataframe, plot_duration, init_time,  plot_start=0.0):

    for signal_name in signal_map:

        signal_raw = dataframe[signal_name][::DECIMATE]
        t = signal_raw.index - init_time
        v = signal_raw.values
        ax.plot(t, v, lw=1.2, c=COL_SIG_1, zorder=4, label=signal_name)

        signal_temp = signal_raw[signal_raw.index > init_time]
        signal_values = signal_temp.values

   
        config_subplot(ax, plot_duration, plot_values=signal_values, signal_name=signal_name, plot_start=plot_start)
        
        
        config_subplot_labels(ax)


    return


def config_subplot(ax, plot_duration, plot_start = 0.0, plot_values = None, signal_name=""):

    # initial_ylims = ax.get_ylim()

    # Grid Properties:
    ax.set_axisbelow(True)
    ax.grid(which='major', linestyle='--', linewidth=0.5, color='black', alpha=0.4, zorder=-100)
    ax.grid(which='minor', linestyle='--', linewidth=0.5, color='black', alpha=0.1, zorder=-100)

    # X-axis
    xlims = [plot_start, plot_duration]

    #  Creates a XMajor Loactor that attempts to place
    xmajor_locator = matplotlib.ticker.MaxNLocator(nbins=5, steps=[1, 2, 4, 5, 10])
    major_xticks_values = xmajor_locator.tick_values(xlims[0], xlims[1])
    ax.set_xticks(major_xticks_values, labels=[])
    ax.xaxis.set_major_locator(xmajor_locator)

    #  Creates a yMajor Loactor that attempts to place
    ymajor_locator = matplotlib.ticker.MaxNLocator(nbins=5, steps=[1, 2, 4, 5, 10])
    ax.yaxis.set_major_locator(ymajor_locator)
    yminor_locator = matplotlib.ticker.MaxNLocator(nbins=20, steps=[1, 2, 4, 5, 10])
    ax.yaxis.set_minor_locator(yminor_locator)

    ax.set_facecolor(SECONDARY_COLOUR)

    # The following minor tick locator is constructed by trying to make less than 20 tickx
    # on the x axis. We try some multiples, 1,2,3,10. And the largest that fits becomes the division.
    minor_tick_options = pd.Series([1, 2, 4, 10])
    rem = minor_tick_options[(minor_tick_options * len(major_xticks_values)) <= 20]
    if len(rem) > 0:
        xminor_locator = matplotlib.ticker.AutoMinorLocator(n=max(rem))
        ax.xaxis.set_minor_locator(xminor_locator)

    #showX-axis labels
    ax.tick_params(axis='x', colors=AXIS_COLOUR, labelsize=8, direction='out', which='both')
    xmajor_formatter = matplotlib.ticker.EngFormatter(unit='s', sep='')
    ax.xaxis.set_major_formatter(xmajor_formatter)

    if plot_values is not None:
        current_y_min = min(plot_values)
        current_y_max = max(plot_values)
    else:
        current_y_min = -1
        current_y_max = 1


    current_y_range = current_y_max - current_y_min
    min_y_range = 0.05
    y_min = current_y_min
    y_max = current_y_max


    # auto scale and add margin
    if signal_name in ["POC_P_MW"]:
        offset = max(current_y_range * 0.2, 0.1)
        y_min = current_y_min - offset
        y_max = current_y_max + offset
        ax.set_ylim(y_min, y_max)
        
    if signal_name in ["BESS1_P_MW", "BESS1_Q_MVAR"]:
        offset = max(current_y_range * 0.2, 0.1)
        y_min = current_y_min - offset
        y_max = current_y_max + offset
        ax.set_ylim(y_min, y_max)
        
    elif signal_name in ["POC_Q_MVAR"]:
        offset = max(current_y_range * 0.2, 0.1)
        y_min = current_y_min - offset
        y_max = current_y_max + offset
        ax.set_ylim(y_min, y_max)
        
    elif signal_name in ["WT1_P_MW", "WT1_Q_MVAR"]:
        offset = max(current_y_range * 0.2, 0.1)
        y_min = current_y_min - offset
        y_max = current_y_max + offset
        ax.set_ylim(y_min, y_max)
        
    elif signal_name in ["POC_V_PU", "WT1_V_PU", "BESS1_V_PU"]:
        offset = max(current_y_range * 0.2, 0.01)
        y_min = current_y_min - offset
        y_max = current_y_max + offset
        ax.set_ylim(y_min, y_max)
        
    elif signal_name in ["POC_ANG_DEG"]:
        offset = max(current_y_range * 0.2, 0.01)
        y_min = current_y_min - offset
        y_max = current_y_max + offset
        ax.set_ylim(y_min, y_max)
        
    elif "_HZ" in signal_name:
        offset = max(min(current_y_range * 0.2, 20), 0.1)
        y_min = current_y_min - offset
        y_max = current_y_max + offset
        ax.set_ylim(y_min, y_max)
        
    else:
        pass

    # Xlims
    ax.set_xlim(xlims[0],xlims[1])
 
    # Y-axis    
    ax.tick_params(axis='y',  colors=AXIS_COLOUR, labelsize=8, direction='in', which='both')

    ax.spines["bottom"].set_color(AXIS_COLOUR)      
    ax.spines["top"].set_color(AXIS_COLOUR)
    ax.spines["left"].set_color(AXIS_COLOUR)
    ax.spines["right"].set_color(AXIS_COLOUR)

def config_subplot_labels(ax,  legend_y=1.12, legend_ncol=3):

    ax.legend(frameon=False, fontsize='x-small', bbox_to_anchor=(1.0, legend_y), ncol=legend_ncol,
            borderpad=0, columnspacing=1, handletextpad=0.3, handlelength=1.2, labelcolor=AXIS_COLOUR)
        
def plot_logo(fig):
    gs_logo = gridspec.GridSpec(6, 6, figure=fig)
    logo_ax = fig.add_subplot(gs_logo[0,-1])
    logo = plt.imread("./resources/windlab.jpg")
    logo_ax.imshow(logo)
    logo_ax.axis("off")
    
    
def save_plot(save_path, png = True, pdf = True):

    path = (str(save_path))
    if pdf:
        plt.savefig(path + ".pdf", bbox_inches='tight', dpi=DPI, format='pdf')
    if png:
        plt.savefig(path + ".png", bbox_inches='tight', dpi=DPI, format='png')
    plt.cla()  # Clear the current axes.
    plt.clf()  # Clear the current figure.
    plt.close()  # Closes all the figure windows.

def mpl_setup(w,h):

    plt.clf()

    plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.labelsize'] = 'medium'
    plt.rcParams['axes.titlesize'] = 'medium'
    plt.rcParams['axes.grid'] = True
    plt.rcParams['axes.grid.axis'] = 'both'
    plt.rcParams['figure.titleweight'] = 'bold'
    plt.rcParams['figure.dpi'] = DPI

    fig = plt.figure(facecolor=PRIMARY_COLOUR) #subplot(number_of_rows, 4)
    
    cm = 1 / 2.54
    fig.set_size_inches(w * cm, h * cm)
    plt.subplots_adjust(left=0.04, right=0.975, bottom=0.05, top=0.96, wspace=0.12, hspace=0.4)


    return fig