import math
from icecream import ic
from rengen.plotting.Plotter import Plotter
import os
import json
import logging
from typing import List, Optional, Tuple, Union, Dict
from pathlib import Path
import pandas as pd
import matplotlib
from rengen.utils.gui_utils import std_script_title
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

LW_NORM = 0.8
LW_REF = 1.2 * LW_NORM

NUM_WTG = 12
NUM_BESS = 2

PLOT_INIT = False

class BBWFPlotter(Plotter):

    def plot_from_df_and_dict(
            self,
            df: pd.DataFrame,
            spec_dict: Dict,
            png_path: Optional[Union[Path, str]] = None,
            pdf_path: Optional[Union[Path, str]] = None,
    ):
        
        
        plt.clf()
        
        if PLOT_INIT:
            self.plot_start = 0
        else:
            self.plot_start = float(spec_dict["substitutions"]["TIME_Full_Init_Time_sec"])
            
        self.plot_duration = df["POC_Vrms_pu"].index[-1] - self.plot_start
        self.plot_end = df["POC_Vrms_pu"].index[-1]
        
        # Balanced or Unbalanced Fault
        unbalanced = False
        if 'Fault_Type_v' in spec_dict:
            try:
                fault_type_numeric_value = pd.to_numeric(spec_dict['Fault_Type_v'])
                fault_type_is_number = np.isfinite(fault_type_numeric_value)
            except (ValueError, TypeError):
                fault_type_is_number = False

            if fault_type_is_number and fault_type_numeric_value != 7:
                unbalanced = True

        # Table Data
        default_table_data = [
            ['Project', "Bungaban Wind Farm"],
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

        # Append Optional Tags to Table Data
        for col_name, label in optional_table_information:
            if col_name in spec_dict:
                if not str(spec_dict[col_name]) == "nan":
                    default_table_data.append([label, f"{spec_dict[col_name]}"])

        table_data = default_table_data
        if "Table_Data" in spec_dict and not pd.isna(spec_dict["Table_Data"]):
            try:
                table_data = json.loads(spec_dict["Table_Data"])
            except json.decoder.JSONDecodeError as e:
                table_data = default_table_data

        # Subplot Layout
        if unbalanced:
            number_of_rows = 9
            number_of_misc_plots = 8
        else:
            number_of_rows = 6
            number_of_misc_plots = 5

        # Page Size and Spacing
        fig = plt.figure(facecolor=PRIMARY_COLOUR)
        fig.set_size_inches(42 / 2.54, 29.7 / 2.54) 
        plt.subplots_adjust(left=0.04, right=0.975, bottom=0.05, top=0.96, wspace=0.12, hspace=0.4)

        # Construct Layout
        layout_main = gridspec.GridSpec(1, 4, figure=fig)
        column_1 = gridspec.GridSpecFromSubplotSpec(number_of_rows, 1, subplot_spec=layout_main[0])
        column_2 = gridspec.GridSpecFromSubplotSpec(number_of_rows, 1, subplot_spec=layout_main[1])
        column_3 = gridspec.GridSpecFromSubplotSpec(number_of_rows, 1, subplot_spec=layout_main[2])
        column_4 = gridspec.GridSpecFromSubplotSpec(number_of_rows, 1, subplot_spec=layout_main[3])
        column_other = gridspec.GridSpecFromSubplotSpec(number_of_misc_plots, 1, subplot_spec=column_4[:number_of_rows - 1])

        # POC Axes
        ax_poc: List[plt.Axes] = []
        for i in range(number_of_rows):
            ax_poc.append(fig.add_subplot(column_1[i, -1]))

        # WTG Axes
        ax_wtg: List[plt.Axes] = []
        for i in range(number_of_rows):
            ax_wtg.append(fig.add_subplot(column_2[i, -1]))

        # BESS Axes
        ax_bess: List[plt.Axes] = []
        for i in range(number_of_rows):
            ax_bess.append(fig.add_subplot(column_3[i, -1]))

        # Table Axis
        ax_table = fig.add_subplot(column_other[0:2, -1])
        ax_table.axis('off')
        
        # inject neg seq current plot for unbalanced fault
        n = 0
        if unbalanced:
            ax_neg_seq_current = fig.add_subplot(column_other[-4, -1])
            n = 1
            
        # ax_33kv_bus = fig.add_subplot(column_other[3+n, -1])
        # ax_oltc = fig.add_subplot(column_other[-4, -1])
        # ax_tap_pos = fig.add_subplot(column_other[-3, -1])
        ax_wtg_frt = fig.add_subplot(column_other[-2,-1])
        ax_bess_frt = fig.add_subplot(column_other[-1,-1])
        ax_frequency = fig.add_subplot(column_4[-1, -1])

        # Column Titles
        ax_poc[0].set_title('POC', style='italic', fontsize='medium', loc='center', y=1.1, pad=5, color=AXIS_COLOUR)
        ax_wtg[0].set_title('WTG', style='italic', fontsize='medium', loc='center', y=1.1, pad=5, color=AXIS_COLOUR)
        ax_bess[0].set_title('BESS', style='italic', fontsize='medium', loc='center', y=1.1, pad=5, color=AXIS_COLOUR)
        ax_table.set_title('Test Parameters', style='italic', fontsize='medium', loc='center', y=1.04, pad=5, color=AXIS_COLOUR)
        


        # Plot Table
        table_colour = np.empty_like(table_data, dtype='object')
        for i, _ in enumerate(table_colour):
            table_colour[i] = [SECONDARY_COLOUR, SECONDARY_COLOUR]

        ax_table.table(cellText=table_data, cellLoc='left', loc='upper right', cellColours=table_colour)

        def set_table_edge_color(ax, color):
            table: Table
            for child in ax.get_children():
                if isinstance(child, Table):
                    table = child
                    break

            cell: Cell
            for cell in table.get_children():
                if isinstance(cell, Cell):
                    cell.set_edgecolor(AXIS_COLOUR)
                    cell.set_linewidth(0.75)
                    cell.get_text().set_color(AXIS_COLOUR)

        if len(table_data) == 0:
            set_table_edge_color(ax_table, AXIS_COLOUR)
            

    

        # PPOC Plotting
        pref_mw = df['POC_WTG_Pref_MW'][self.plot_start:self.plot_end][::DECIMATE]
        ppoc_mw = df['POC_P_MW'][self.plot_start:self.plot_end][::DECIMATE]
        poc_fdroop_mw = df['POC_Pref_Fdroop_MW'][self.plot_start:self.plot_end][::DECIMATE]
        self.signal_plot(
            ax=ax_poc[0],
            title='POC: P [MW]',
            traces=[
                    ('Pref + Fdroop', LW_NORM, (0.1,0.7,0.1,0.6), 0, poc_fdroop_mw), 
                    ('Pref', LW_REF, COL_REF, 1, pref_mw),
                    ('Ppoc', LW_NORM, COL_POC, 2, ppoc_mw), ],
            min_y_range=40,
        )

        # QPOC Plotting
        qpoc_mvar = df['POC_Q_MVAr'][self.plot_start:self.plot_end][::DECIMATE]
        qref_mvar = df['Qref_droop_MVAr'][self.plot_start:self.plot_end][::DECIMATE]
        ideal_qref_mvar = df['Qref_droop_Ideal_MVAr'][self.plot_start:self.plot_end][::DECIMATE]
        self.signal_plot(
            ax=ax_poc[1],
            title='POC: Q [MVAr]',
            traces=[('Qref + Vdroop', LW_REF, (0.1,0.7,0.1,0.6), 1, ideal_qref_mvar),
                    ('Qref', LW_REF, COL_REF, 1, qref_mvar),
                    ('Qpoc', LW_NORM, COL_POC, 2, qpoc_mvar), ],
            min_y_range=20,
        )

        # VPOC Plotting
        vpoc_pu = df['POC_Vrms_pu'][self.plot_start:self.plot_end][::DECIMATE]
        vref_pu = df['POC_Vref_pu'][self.plot_start:self.plot_end][::DECIMATE]
        self.signal_plot(
            ax=ax_poc[2],
            title='POC: V [.pu]',
            traces=[('Vref', LW_REF, COL_REF, 1, vref_pu),
                    ('Vpoc', LW_NORM, COL_POC, 2, vpoc_pu), ],
            min_y_range=0.05
        )

        if unbalanced:
            try:
                poc_i_neg_pu = df['POC_I_neg_pu'][self.plot_start:self.plot_end][::DECIMATE]
                self.signal_plot(
                    ax=ax_neg_seq_current,
                    title='POC Neg. Seq. Current: pu',
                    traces=[('POC', LW_NORM, COL_POC, 1, poc_i_neg_pu)],
                    min_y_range=0.1,
                )
            except KeyError:
                print("POC_I_neg_pu not found in psout/pkl.")
                pass


            # VPOC ABC Phase Voltages
            poc_vrms_a = df['POC_Va_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]
            poc_vrms_b = df['POC_Vb_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]
            poc_vrms_c = df['POC_Vc_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]

            ymins = []
            ymaxs = []

            for row, phase, signal in zip(
                    [3, 4, 5],
                    ['A', 'B', 'C'],
                    [poc_vrms_a, poc_vrms_b, poc_vrms_c]):
                ax = ax_poc[row]

                self.signal_plot(
                    ax=ax,
                    title='POC: Phase ' + phase + ' Voltage [.pu]',
                    traces=[('POC', LW_NORM, COL_POC, 1, signal)],
                    min_y_range=0.1,
                )

                ymin, ymax = ax.get_ylim()
                ymins.append(ymin)
                ymaxs.append(ymax)

            for row in [3, 4, 5]:
                ax = ax_poc[row]
                ax.set_ylim(min(ymins), max(ymaxs))

        # POC ID Plotting
        id_poc_pu = df['POC_Id_pos_pu'][self.plot_start:self.plot_end][::DECIMATE]
        self.signal_plot(
            ax=ax_poc[number_of_rows - 3],
            title='POC: Id [.pu]',
            traces=[('POC', LW_NORM, COL_POC, 1, id_poc_pu)],
            min_y_range=0.05,
        )

        # POC Iq Plotting
        iq_poc_pu = df['POC_Iq_pos_pu'][self.plot_start:self.plot_end][::DECIMATE]
        self.signal_plot(
            ax=ax_poc[number_of_rows - 2],
            title='POC: Iq [.pu]',
            traces=[('POC', LW_NORM, COL_POC, 1, iq_poc_pu)],
            min_y_range=0.05,
        )
        
        # POC Angle Plotting
        self.signal_plot(
            ax=ax_poc[-1],
            title='POC: Angle [deg]',
            traces=[("POC", LW_NORM, COL_POC, 5, df['POC_Angle_deg'][self.plot_start:self.plot_end][::DECIMATE]), ],
            min_y_range=3,
            time_axis_on=True,
        )
    
        # POC Frequency Plotting
        self.signal_plot(
            ax=ax_frequency,
            title='POC: Frequency [Hz]',
            traces=[("POC", LW_NORM, COL_POC, 5, df['POC_Freq_Hz'][self.plot_start:self.plot_end][::DECIMATE]), ],
            min_y_range=0.5,
            time_axis_on=True,
        )
        
            
        # FRT and Trip plots
        traces = []
        for i in range(1,NUM_WTG+1):
            traces.append(("", LW_NORM, COL_LVRT, i, df[f'WTG_LVRT_Flag:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
            traces.append(("", LW_NORM, COL_HVRT, i, df[f'WTG_HVRT_Flag:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
            traces.append(("", LW_NORM, COL_TRIP, i, df[f'WTG_Trip_Flag:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
            
        ax_wtg_frt.plot([], [], c=COL_LVRT, label=f"LVRT WTG1-{NUM_WTG}")
        ax_wtg_frt.plot([], [], c=COL_HVRT, label=f"HVRT WTG1-{NUM_WTG}")
        ax_wtg_frt.plot([], [], c=COL_TRIP, label=f"TRIP WTG1-{NUM_WTG}")
        
        self.signal_plot(
            ax=ax_wtg_frt,
            title='WTG FRT and Trip Flags',
            traces=traces,
            min_y_range=3,        )
        
        traces = []
        for i in range(1,NUM_BESS+1):
            traces.append(("", LW_NORM, COL_LVRT, i, df[f'BESS_LVRT_Flag:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
            traces.append(("", LW_NORM, COL_HVRT, i, df[f'BESS_HVRT_Flag:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
            traces.append(("", LW_NORM, COL_TRIP, i, df[f'BESS_Trip_Flag:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
            
        ax_bess_frt.plot([], [], c=COL_LVRT, label=f"LVRT BESS1-{NUM_BESS}")
        ax_bess_frt.plot([], [], c=COL_HVRT, label=f"HVRT BESS1-{NUM_BESS}")
        ax_bess_frt.plot([], [], c=COL_TRIP, label=f"TRIP BESS1-{NUM_BESS}")
        
        self.signal_plot(
            ax=ax_bess_frt,
            title='BESS FRT and Trip Flags',
            traces=traces,
            min_y_range=3,
        )
        

        # WTG P Plotting
        wtg_p_traces=[]
        wtg_p_traces.append(('Pref', LW_REF, COL_REF, 20, df['WTG_Unit_Pref_MW'][self.plot_start:self.plot_end][::DECIMATE]),)
        for i in range(1,NUM_WTG+1):
            wtg_p_traces.append(("", LW_NORM, COL_SIG_2, i, df[f'WTG_Terminal_P_MW:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
        ax_wtg[0].plot([], [], c=COL_SIG_2, label=f"WTG1-{NUM_WTG}")
        self.signal_plot(
            ax=ax_wtg[0],
            title='WTG: P [MW]',
            traces=wtg_p_traces,
            min_y_range=1,
        )
        

        # WTG Q Plotting
        wtg_q_traces=[]
        # wtg_q_traces.append(('Qref', LW_REF, COL_REF, 20, df['WTG_Unit_Qset_MVAr:1'][self.plot_start:self.plot_end][::DECIMATE]),)
        for i in range(1,NUM_WTG+1):
            wtg_q_traces.append(("", LW_NORM, COL_SIG_2, i, df[f'WTG_Terminal_Q_MVAr:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
        ax_wtg[1].plot([], [], c=COL_SIG_2, label=f"WTG1-{NUM_WTG}")
        self.signal_plot(
            ax=ax_wtg[1],
            title='WTG: Q [MVAr]',
            traces=wtg_q_traces,
            min_y_range=1,
        )

        # WTG V Plotting
        wtg_v_traces=[]
        wtg_v_traces.append(('Vref', LW_REF, COL_REF, 20, df['WTG_Vset_pu:1'][self.plot_start:self.plot_end][::DECIMATE]),)
        for i in range(1,NUM_WTG+1):
            wtg_v_traces.append(("", LW_NORM, COL_SIG_2, i, df[f'WTG_Terminal_V_pu:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
        ax_wtg[2].plot([], [], c=COL_SIG_2, label=f"WTG1-{NUM_WTG}")
        self.signal_plot(
            ax=ax_wtg[2],
            title='WTG Terminal: V [.pu]',
            traces=wtg_v_traces,
            min_y_range=0.05,
        )

        if unbalanced:
            # WTG V ABC PHASE Plotting
            traces_a = []
            traces_b = []
            traces_c = []
            for i in range(1, NUM_WTG + 1):
                # signal_label, lw, colour, order, signal
                traces_a.append(("", LW_NORM, COL_SIG_2, i,
                                     df[f'WTG_Va_rms_pu:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
                traces_b.append(("", LW_NORM, COL_SIG_2, i,
                                     df[f'WTG_Vb_rms_pu:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
                traces_c.append(("", LW_NORM, COL_SIG_2, i,
                                     df[f'WTG_Vc_rms_pu:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
                ymins = []
                ymaxs = []
                for row, phase, traces in zip(
                        [3, 4, 5],
                        ['A', 'B', 'C'],
                        [traces_a, traces_b, traces_c]):
                    ax = ax_wtg[row]

                    self.signal_plot(
                        ax=ax,
                        title='WTG: Phase ' + phase + ' Voltage [.pu]',
                        traces=traces,
                        min_y_range=0.1,
                    )

                    ymin, ymax = ax.get_ylim()
                    ymins.append(ymin)
                    ymaxs.append(ymax)

                for row in [3, 4, 5]:
                    ax = ax_wtg[row]
                    # ax.hlines(lvrt_th_in, 0, plot_duration, linestyle='-', lw=0.5, color='0.5', zorder=1, label="FRT-in")
                    # ax.hlines(hvrt_th_in, 0, plot_duration, linestyle='-', lw=0.5, color='0.5', zorder=1)
                    # ax.hlines(lvrt_th_out, 0, plot_duration, linestyle='-.', lw=0.5, color='0.5', zorder=1)
                    # ax.hlines(hvrt_th_out, 0, plot_duration, linestyle='-.', lw=0.5, color='0.5', zorder=1, label="FRT-out")
                    ax.set_ylim(min(ymins), max(ymaxs))

        # WTG Id Plotting
        wtg_id_traces=[]
        for i in range(1,NUM_WTG+1):
            wtg_id_traces.append(("", LW_NORM, COL_SIG_2, i, df[f'WTG_Id_pu:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
        ax_wtg[number_of_rows - 3].plot([], [], c=COL_SIG_2, label=f"WTG1-{NUM_WTG}")
        self.signal_plot(
            ax=ax_wtg[number_of_rows - 3],
            title='WTG: Id [.pu]',
            traces=wtg_id_traces,
            min_y_range=0.05,
        )

        # WTG Iq Plotting
        wtg_iq_traces=[]
        for i in range(1,NUM_WTG+1):
            wtg_iq_traces.append(("", LW_NORM, COL_SIG_2, i, df[f'WTG_Iq_pu:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
        ax_wtg[number_of_rows - 2].plot([], [], c=COL_SIG_2, label=f"WTG1-{NUM_WTG}")
        self.signal_plot(
            ax=ax_wtg[number_of_rows - 2],
            title='WTG: Iq [.pu]',
            traces=wtg_iq_traces,
            min_y_range=0.05,
        )
        
        # WTG Angle Plotting
        wtg_angle_traces=[]
        for i in range(1,NUM_BESS+1):
            wtg_angle_traces.append(("", LW_NORM, COL_SIG_2, i, df[f'WTG_Angle_deg:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
        ax_wtg[-1].plot([], [], c=COL_SIG_2, label=f"WTG1-{NUM_WTG}")
        self.signal_plot(
            ax=ax_wtg[-1],
            title='WTG: Angle [deg]',
            traces=wtg_angle_traces,
            min_y_range=3,
            time_axis_on=True,
        )

        # BESS P Plotting
        bess_p_traces=[]
        bess_p_traces.append(('Pref', LW_REF, COL_REF, 20, df['BESS_Unit_Pref_MW:1'][self.plot_start:self.plot_end][::DECIMATE]),)
        for i in range(1,NUM_BESS+1):
            bess_p_traces.append(("", LW_NORM, COL_SIG_2, i, df[f'BESS_Terminal_P_MW:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
        ax_bess[0].plot([], [], c=COL_SIG_2, label=f"BESS1-{NUM_BESS}")
        self.signal_plot(
            ax=ax_bess[0],
            title='BESS: P [MW]',
            traces=bess_p_traces,
            min_y_range=1,
        )

        # BESS Q Plotting
        bess_q_traces=[]
        bess_q_traces.append(('Qref', LW_REF, COL_REF, 20, df['BESS_Unit_Qset_MVAr:1'][self.plot_start:self.plot_end][::DECIMATE]),)
        for i in range(1,NUM_BESS+1):
            bess_q_traces.append(("", LW_NORM, COL_SIG_2, i, df[f'BESS_Terminal_Q_MVAr:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
        ax_bess[1].plot([], [], c=COL_SIG_2, label=f"BESS1-{NUM_BESS}")
        self.signal_plot(
            ax=ax_bess[1],
            title='BESS: Q [MVAr]',
            traces=bess_q_traces,
            min_y_range=1,
        )

        # BESS V Plotting
        bess_v_traces=[]
        bess_v_traces.append(('Vref', LW_REF, COL_REF, 20, df['BESS_Vset_pu:1'][self.plot_start:self.plot_end][::DECIMATE]),)
        for i in range(1,NUM_BESS+1):
            bess_v_traces.append(("", LW_NORM, COL_SIG_2, i, df[f'BESS_Terminal_V_pu:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
        ax_bess[2].plot([], [], c=COL_SIG_2, label=f"BESS1-{NUM_BESS}")
        self.signal_plot(
            ax=ax_bess[2],
            title='BESS Terminal: V [.pu]',
            traces=bess_v_traces,
            min_y_range=0.05,
        )

        # # BESS V ABC PHASE Plotting
        # if unbalanced:

        #     traces_a = [('B1', LW_NORM, COL_SIG_2, 4, df['BESS1_Va_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]),
        #                 ('B2', LW_NORM, COL_SIG_2, 3, df['BESS2_Va_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]),
        #                 ('B3', LW_NORM, COL_SIG_3, 2, df['BESS3_Va_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]),
        #                 ('B4', LW_NORM, COL_SIG_4, 1, df['BESS4_Va_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]), ]

        #     traces_b = [('B1', LW_NORM, COL_SIG_2, 4, df['BESS1_Vb_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]),
        #                 ('B2', LW_NORM, COL_SIG_2, 3, df['BESS2_Vb_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]),
        #                 ('B3', LW_NORM, COL_SIG_3, 2, df['BESS3_Vb_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]),
        #                 ('B4', LW_NORM, COL_SIG_4, 1, df['BESS4_Vb_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]), ]

        #     traces_c = [('B1', LW_NORM, COL_SIG_2, 4, df['BESS1_Vc_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]),
        #                 ('B2', LW_NORM, COL_SIG_2, 3, df['BESS2_Vc_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]),
        #                 ('B3', LW_NORM, COL_SIG_3, 2, df['BESS3_Vc_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]),
        #                 ('B4', LW_NORM, COL_SIG_4, 1, df['BESS4_Vc_rms_pu'][self.plot_start:self.plot_end][::DECIMATE]), ]

        #     ymins = []
        #     ymaxs = []
        #     for row, phase, traces in zip(
        #             [3, 4, 5],
        #             ['A', 'B', 'C'],
        #             [traces_a, traces_b, traces_c]):
        #         ax = ax_bess[row]

        #         signal_plot(
        #             ax=ax,
        #             title='BESS: Phase ' + phase + ' Voltage [.pu]',
        #             traces=traces,
        #             start=plot_start,
        #             duration=plot_duration,
        #             min_y_range=0.1,
        #         )

        #         ymin, ymax = ax.get_ylim()
        #         ymins.append(ymin)
        #         ymaxs.append(ymax)

        #     for row in [3, 4, 5]:
        #         ax = ax_bess[row]
        #         ax.set_ylim(min(ymins), max(ymaxs))

        # BESS Id Plotting
        bess_id_traces=[]
        for i in range(1,NUM_BESS+1):
            bess_id_traces.append(("", LW_NORM, COL_SIG_2, i, df[f'BESS_Id_pu:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
        ax_bess[number_of_rows - 3].plot([], [], c=COL_SIG_2, label=f"BESS1-{NUM_BESS}")
        self.signal_plot(
            ax=ax_bess[number_of_rows - 3],
            title='BESS: Id [.pu]',
            traces=bess_id_traces,
            min_y_range=0.05,
        )

        # BESS Iq Plotting
        bess_iq_traces=[]
        for i in range(1,NUM_BESS+1):
            bess_iq_traces.append(("", LW_NORM, COL_SIG_2, i, df[f'BESS_Iq_pu:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
        ax_bess[number_of_rows - 2].plot([], [], c=COL_SIG_2, label=f"BESS1-{NUM_BESS}")
        self.signal_plot(
            ax=ax_bess[number_of_rows - 2],
            title='BESS: Iq [.pu]',
            traces=bess_iq_traces,
            min_y_range=0.05,
        )

        
        # BESS Angle Plotting
        bess_angle_traces=[]
        for i in range(1,NUM_BESS+1):
            bess_angle_traces.append(("", LW_NORM, COL_SIG_2, i, df[f'BESS_Angle_deg:{i}'][self.plot_start:self.plot_end][::DECIMATE]))
        ax_bess[-1].plot([], [], c=COL_SIG_2, label=f"BESS1-{NUM_BESS}")
        self.signal_plot(
            ax=ax_bess[-1],
            title='BESS: Angle [deg]',
            traces=bess_angle_traces,
            min_y_range=3,
            time_axis_on=True,
        )       
        
        # -----------------------------------------------------------------------------------------------------------------
        #  Analysis Annotations: 
        # -----------------------------------------------------------------------------------------------------------------
        
        if "expected_poc_p_ramp_t" in spec_dict["analysis"]:
            
            self.plot_curve(
                x_values = [float(v) for v in list(spec_dict["analysis"]["expected_poc_p_ramp_t"])],
                y_values =  [float(v) for v in list(spec_dict["analysis"]["expected_poc_p_ramp_mw"])],
                ax=ax_poc[0],
                label="Expected POC P NEMDE Ramp"
            )
        
        if "Fault_iq_AEMO_Rise_Time" in spec_dict["analysis"]:
            
            rise_time = float(spec_dict["analysis"]["Fault_iq_AEMO_Rise_Time"])
            self.plot_text_annotations(
                text=f"Rise Time = {1000 * rise_time:.2f} ms", 
                x_pos=rise_time+float(spec_dict["Fault_Time"]), 
                y_pos_ratio=0.1, 
                ax=ax_poc[number_of_rows - 2]
            )
            self.plot_vertical_lines(
                x_pos=float(spec_dict["Fault_Time"]), 
                ax=ax_poc[number_of_rows - 2]
            )
            self.plot_vertical_lines(
                x_pos=float(spec_dict["analysis"]["Fault_iq_AEMO_Rise_Time"])+float(spec_dict["Fault_Time"]), 
                ax=ax_poc[number_of_rows - 2]
            )
            self.plot_horrisontal_lines(
                y_pos=float(spec_dict["analysis"]["Fault_iq_mcc_pu"]), 
                ax=ax_poc[number_of_rows - 2]
            )
            
            
        if "Fault_iq_AEMO_Settling_Time" in spec_dict["analysis"]:   
            
            settling_time = float(spec_dict["analysis"]["Fault_iq_AEMO_Settling_Time"])
            self.plot_text_annotations(
                text=f"Settling Time = {1000 * settling_time:.2f} ms", 
                x_pos=settling_time+float(spec_dict["Fault_Time"]), 
                y_pos_ratio=0.9, 
                ax=ax_poc[number_of_rows - 2]
            )
            self.plot_vertical_lines(
                x_pos=float(spec_dict["analysis"]["Fault_iq_AEMO_Settling_Time"])+float(spec_dict["Fault_Time"]), 
                ax=ax_poc[number_of_rows - 2])         
            
            
        if "active_power_rampdown_start_s" in spec_dict["analysis"]: 
            
            self.plot_text_annotations(
                text="POC P Ramp Down Region", 
                x_pos=float(spec_dict["analysis"]["active_power_rampdown_start_s"]),
                y_pos_ratio=0.1, 
                ax=ax_poc[0]
            )
            self.plot_vertical_lines(
                x_pos=float(spec_dict["analysis"]["active_power_rampdown_start_s"]), 
                ax=ax_poc[0]
            )
            self.plot_vertical_lines(
                x_pos=float(spec_dict["analysis"]["active_power_rampdown_expected_end_s"]), 
                ax=ax_poc[0]
            )
            self.plot_horrisontal_lines(
                y_pos=float(spec_dict["analysis"]["active_power_rampdown_threshold_mw"]), 
                ax=ax_poc[0]
            )
        
        
        
        # -----------------------------------------------------------------------------------------------------------------
        #  Save Plot
        # -----------------------------------------------------------------------------------------------------------------

        fig.align_labels()
        
        plt.savefig(pdf_path, bbox_inches='tight', dpi=300, format='pdf')
        plt.savefig(png_path, bbox_inches='tight', dpi=300, format='png')
        
        plt.cla() 
        plt.clf()
        plt.close()
        
        return

    def plot_curve(self, x_values: List[float], y_values: List[float], ax: plt.Axes, label=None):

        # ymin, ymax = ax.get_ylim()
        # xmin, xmax = ax.get_xlim()
        ax.plot(x_values, y_values, color=(0.3,0.3,0.3,0.6), linestyle='--', lw=LW_REF, label=label)
        ax.legend(frameon=False, fontsize='x-small', bbox_to_anchor=(1.0, 0.0), ncol=6,
                    borderpad=0, columnspacing=1, handletextpad=0.3, handlelength=1.2, labelcolor=AXIS_COLOUR)
        # ax.set_ylim(ymin, ymax)
        # ax.set_xlim(xmin, xmax)

    def plot_vertical_lines(self, x_pos: float, ax: plt.Axes): 
        ymin, ymax = ax.get_ylim()
        ax.vlines(x_pos,ymin,ymax,colors=(0.3,0.3,0.3,0.6),linestyle='--',lw=LW_REF)
        ax.legend(frameon=False, fontsize='x-small', bbox_to_anchor=(1.0, 0.0), ncol=6,
                    borderpad=0, columnspacing=1, handletextpad=0.3, handlelength=1.2, labelcolor=AXIS_COLOUR)
    
    
    def plot_horrisontal_lines(self, y_pos: float, ax: plt.Axes): 
        xmin, xmax = ax.get_xlim()
        ax.hlines(y_pos,xmin,xmax,colors=(0.3,0.3,0.3,0.6),linestyle='--',lw=LW_REF)
        ax.legend(frameon=False, fontsize='x-small', bbox_to_anchor=(1.0, 0.0), ncol=6,
                    borderpad=0, columnspacing=1, handletextpad=0.3, handlelength=1.2, labelcolor=AXIS_COLOUR)
    
    def plot_text_annotations(self, text: str, x_pos: float, y_pos_ratio: float, ax: plt.Axes):
        ymin, ymax = ax.get_ylim()
        xmin, xmax = ax.get_xlim()
        x_pad = (xmax-xmin)*0.02
        y_pos = ymin + y_pos_ratio*(ymax-ymin)
        ax.text(x_pos+x_pad, y_pos, text, fontsize='x-small', color=AXIS_COLOUR)    
        


    def configure_subplot(self, ax, subplot_title, legendon=True, xlabelson=False, min_y_range=None,
                        fixed_y_lims=None):
        # Grid Properties:
        ax.set_axisbelow(True)
        ax.grid(which='major', linestyle='--', linewidth=0.5, color='darkgray', alpha=0.8, zorder=-100)
        ax.grid(which='minor', linestyle='--', linewidth=0.5, color='darkgray', alpha=0.4, zorder=-100)

        # X-Axis Propertie:
        ax.set_xlim(0, self.plot_duration)
        #  Creates a XMajor Loactor that attempts to place
        xmajor_locator = matplotlib.ticker.MaxNLocator(nbins=5, steps=[1, 2, 4, 5, 10])
        major_xticks_values = xmajor_locator.tick_values(0, self.plot_duration)
        ax.set_xticks(major_xticks_values, labels=[])
        ax.xaxis.set_major_locator(xmajor_locator)

        ax.set_facecolor(SECONDARY_COLOUR)

        # The following minor tick locator is constructed by trying to make less than 20 tickx
        # on the x axis. We try some multiples, 1,2,3,10. And the largest that fits becomes the division.
        minor_tick_options = pd.Series([1, 2, 4, 10])
        rem = minor_tick_options[(minor_tick_options * len(major_xticks_values)) <= 20]
        if len(rem) > 0:
            xminor_locator = matplotlib.ticker.AutoMinorLocator(n=max(rem))
            ax.xaxis.set_minor_locator(xminor_locator)

        if xlabelson:
            ax.tick_params(axis='x', colors=AXIS_COLOUR, labelsize=8, direction='out', which='both')
            xmajor_formatter = matplotlib.ticker.EngFormatter(unit='s', sep='')
            ax.xaxis.set_major_formatter(xmajor_formatter)
        else:
            # ax.xaxis.set_major_formatter(matplotlib.ticker.NullFormatter)
            ax.tick_params(axis='x', colors=AXIS_COLOUR, direction='in', which='both')

        ax.tick_params(axis='y', colors=AXIS_COLOUR, labelsize=7, direction='in', which='both')

        # Y - Axis.
        # ymajor_formatter = matplotlib.ticker.EngFormatter(unit='', sep='')
        # ax.yaxis.set_major_formatter(ymajor_formatter)

        if fixed_y_lims is None:
            ax.autoscale(enable=True, axis='y')
            ymin, ymax = ax.get_ylim()
            ax.set_ylim(ymin, ymax)

            if min_y_range is not None:
                ymin, ymax = ax.get_ylim()
                if (ymax - ymin) <= min_y_range:
                    ax.autoscale(enable=True, axis='y', tight=True)
                    ymin, ymax = ax.get_ylim()
                    mid_point = ymin + (ymax - ymin) / 2
                    new_ymax = mid_point + (min_y_range / 2)
                    new_ymin = mid_point - (min_y_range / 2)
                    ax.set_ylim(new_ymin, new_ymax)
        else:
            ax.set_ylim(*fixed_y_lims)
        ax.locator_params(axis='y', nbins=5)
        # Title
        ymin, ymax = ax.get_ylim()
        # ax.text(0.0, ymax, subplot_title, fontsize='small', va='bottom', fontweight='normal')
        ax.set_title(subplot_title, fontsize='small', loc='left', fontweight='normal', y=1, pad=3, color=AXIS_COLOUR)

        # Only ad legend if there are labels.
        _, tmplables = ax.get_legend_handles_labels()
        if tmplables and legendon:
            ax.legend(frameon=False, fontsize='x-small', bbox_to_anchor=(1.0, 0.0), ncol=6,
                    borderpad=0, columnspacing=1, handletextpad=0.3, handlelength=1.2, labelcolor=AXIS_COLOUR)

        ax.spines["bottom"].set_color(AXIS_COLOUR)
        ax.spines["top"].set_color(AXIS_COLOUR)
        ax.spines["left"].set_color(AXIS_COLOUR)
        ax.spines["right"].set_color(AXIS_COLOUR)


    def signal_plot(self, ax, title, traces, fixed_y_lims=None, min_y_range=None, time_axis_on=False):
        
        
        if time_axis_on:
            legendon = False
            xlabelson = True
        else:
            legendon = True
            xlabelson = False

        for signal_label, lw, colour, order, signal in traces:
            if signal_label == "":
                ax.plot(signal.index - self.plot_start, signal.values,
                        lw=lw, c=colour, zorder=order)
            else:
                ax.plot(signal.index - self.plot_start, signal.values,
                        lw=lw, c=colour, zorder=order, label=signal_label)

        if fixed_y_lims is not None:
            self.configure_subplot(ax, title, fixed_y_lims=fixed_y_lims, legendon=legendon, xlabelson=xlabelson)
        elif min_y_range is not None:
            self.configure_subplot(ax, title, min_y_range=min_y_range, legendon=legendon, xlabelson=xlabelson)
        else:
            self.configure_subplot(ax, title, legendon=legendon, xlabelson=xlabelson)




    def read_flag_signal(self, flag_signal):
        """
        Description:

        Function to read and interpret flag signals, returning
        a list of intervals at which the flag is high or low

        Inputs: Flag signal that is either 0, 1 or -1

        Outputs: returns two lists containing the interval when the flag is low and high.
        """

        low_flag = []
        high_flag = []

        previous_value = 0.0
        low_flag_start = 0
        low_flag_stop = 0
        high_flag_start = 0
        high_flag_stop = 0

        for [value, time] in zip(np.append(flag_signal.values, 0.0),
                                np.append(flag_signal.index, flag_signal.index[-1])):

            difference = value - previous_value

            if difference < 0 and value == -1:
                low_flag_start = time

            if difference > 0 and previous_value == -1:
                low_flag_stop = time
                low_flag.append(pd.Series(data=[1.0, 1.0], index=[low_flag_start, low_flag_stop]))

            if difference > 0 and value == 1:
                high_flag_start = time

            if difference < 0 and previous_value == 1:
                high_flag_stop = time
                high_flag.append(pd.Series(data=[1.0, 1.0], index=[high_flag_start, high_flag_stop]))

            previous_value = value

        return low_flag, high_flag


    def plot_flag(self, axs, plot_start, flag_signal, vertical_pos, flag_color):
        """
        Description:

        Function used to plot flags on a given axis
        """

        axs.plot(
            flag_signal.index - plot_start,
            0 * flag_signal.values + vertical_pos,
            linestyle=None,
            marker='|',
            markeredgewidth=0.6,
            markersize=6,
            lw=4,
            solid_capstyle='butt',
            c=flag_color,
            zorder=2,
        )


    def filtered_y_limits(self, signal, ref, padding=0.2, min_y_range=40):
        """
        Description:

        This function is used to return y-axis limits for a plot that contains both a
        measured signal and a reference signal.
        This function Filters the reference signal to ignore instantaneous spikes and
        retain steady changes in reference including saturation.
        """

        ref_filtered = self.low_pass_filter(signal=ref, cut_off=1).values

        y_min = min(min(ref_filtered), min(signal))
        y_max = max(max(ref_filtered), max(signal))

        y_range = y_max - y_min

        y_lim_min = y_min - padding * y_range
        y_lim_max = y_max + padding * y_range

        y_range = y_lim_max - y_lim_min

        if y_range <= min_y_range:
            mid_point = y_lim_min + y_range / 2
            y_lim_max = mid_point + (min_y_range / 2)
            y_lim_min = mid_point - (min_y_range / 2)

        return y_lim_min, y_lim_max


    def low_pass_filter(self, signal: pd.Series, cut_off: float):
        fs = 1 / (signal.index[1] - signal.index[0])
        b, a = butter(1, cut_off, fs=fs, btype='low', analog=False)
        filtered_signal = lfilter(b, a, signal.values - signal.values[0]) + signal.values[0]

        return pd.Series(
            data=filtered_signal,
            index=signal.index,
        )


    