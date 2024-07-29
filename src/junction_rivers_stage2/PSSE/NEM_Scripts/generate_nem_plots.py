import os
import warnings

import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.table import Cell, Table
from scipy.signal import butter, lfilter
# from nem_analysis import nem_analysis, nem_analysis_s52513
import matplotlib.ticker as ticker

DEBUG = False
VERBOSE = True

OUTPUT_DIR = r"G:\Bungaban\PSSE_NEM_Models\Bungaban_PSSE_NEM_Model_v0.1"
CSV_FILE_PATH=r"G:\Bungaban\PSSE_NEM_Models\Bungaban_PSSE_NEM_Model_v0.1\Test_Contingency\studies__psse_nem_network_capabilities.csv"
# S52513_CSV = r"studies__psse_nem_network_capabilities\studies__NEM_s52513.csv"

TITLE = "NEM Study - Case 1"

RESULT_DIRS = [
    ("Test_Contingency","Test_Contingency"),
    # ("LL10_OOS","LL10_Pmax_IS_Night"),
    # ("LL10_Pmax_OOS_Night","LL10_Pmax_IS_Night"),
    # ("LL10_OOS","LL10_Pmin_IS_Day"),
    # ("HL90_OOS","HL90_Pmax_IS_Day"),
    # ("HL90_OOS","HL90_Pmax_IS_Night"),
    # ("HL90_Pmax_OOS_Night","HL90_Pmax_IS_Night"),
    # ("HL90_OOS","HL90_Pmin_IS_Day"),
]

# S5255 Signals

# POC_SIGNALS = [
#     ["POC_P_MW"],
#     ["POC_Q_MVAR"],
#     ["POC_V_PU"],
# ]

# TERMINAL_SIGNALS = [
#     ["WT1_P_MW","BESS1_P_MW"],
#     ["WT1_Q_MVAR","BESS1_Q_MVAR"],
#     ["WT1_V_PU","BESS1_V_PU"],
# ]
# OTHER_SIGNALS = [
#     ["FPOC_HZ"],
#     ["WT1_ID_PU","BESS1_ID_PU"],
#     ["WT1_IQ_PU","BESS1_IQ_PU"],
# ]

# S52512 Signals

ACTIVE_POWER_SIGNALS = [
    ["WSTDN_CLMB_P_MW"],
    ["HALYS_WSTDN_P_MW"],
    ["WSTDN_CLMB_Q_MVAR"],
]
SYNC_MACHINE_ANGLE_SIGNALS = [
    ["KOGAN_ANGLE_DEG"],
    ["STANWELL_ANGLE_DEG"],
    ["TARONG_NORTH_ANGLE_DEG"],
]
BUS_VOLTAGE_SIGNALS_1 = [
    ["WDSTH_275KV_V_PU"],
    ["CLMB_275KV_V_PU"],
    ["WSTDN_275KV_V_PU"],
]
BUS_VOLTAGE_SIGNALS_2 = [
    ["HALYS_275KV_V_PU"],
    ["TARNG_275KV_V_PU"],
]





# Plotting Parameters
PSSE_DECIMATE = 10
DPI = 150

SIGNAL_LINE_WIDTH = 0.8

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

SIGNAL_COLS = [
    COL_SIG_1,
    COL_SIG_2,
    COL_SIG_3,
    COL_SIG_4,
]

DECIMATE = 1
class NEMPostprocessUtility:

    def __init__(self) -> None:
        self.run_nem_plot()
        return


    def run_nem_plot(self, csv_file_path=CSV_FILE_PATH, output_dir=OUTPUT_DIR):

        s52512_csv_data = pd.read_csv(csv_file_path)
        # s52513_csv_data = pd.read_csv(S52513_CSV)

        # s52512_analysis_df = s52512_csv_data.copy()
        # s52513_analysis_df = s52513_csv_data.copy()


        # # Loop through S52513 Results
        # for index, scenario in s52513_csv_data.iloc[::-1].iterrows():

        #     # Get file name from scenario
        #     filename = scenario["File_Name"]

        #     if VERBOSE: print(f"\nFinding results for: {filename}")

        #     for _, is_dir in RESULT_DIRS:

        #             #making plots directory
        #             plot_dir = os.path.join(output_dir, is_dir, "plots")
        #             if not os.path.exists(plot_dir):
        #                 os.makedirs(plot_dir)

        #             # Directories containing pickle files
        #             is_pkl_dir = os.path.join(output_dir, is_dir, "pickle")

        #             # Locate pickle file that matches filename
        #             nem_is_pkl_file = self.find_pickle(is_pkl_dir, filename)

        #             # Skip if no pickle file found
        #             if nem_is_pkl_file == None:
        #                 continue

        #             if VERBOSE: print("\n", OUTPUT_DIR, is_dir, "plots", nem_is_pkl_file[:-4] + "_test_plot")

        #             # Load in service and out of service data
        #             nem_is_df = pd.read_pickle(os.path.join(is_pkl_dir,nem_is_pkl_file))

        #             test_id = nem_is_pkl_file.replace("PSSE_s5.2.5.12","")[:-4]
        #             load_case = is_dir

        #             s52513_file_name = f"{filename}.pdf"

        #             s52513_save_path = os.path.join(output_dir, is_dir, "plots", s52513_file_name)

        #             # Perform Analysis on Scenario
        #             nem_analysis_s52513(
        #                 scenario=scenario,
        #                 results_df=nem_is_df,
        #             )


        #             # Add Columns to Results if analysis produces new columns.
        #             for col in scenario.index:
        #                 if col not in s52513_analysis_df.columns:
        #                     s52513_analysis_df[col] = ""

        #             # Save Results.
        #             s52513_analysis_df.loc[index] = scenario

        #             # Generate Plots
        #             if VERBOSE: print("\nGenerating S52513 Plot...")
        #             self.s52513_plot("S52513 NEM Plot", load_case, test_id, scenario,nem_is_df, s52513_save_path)


        # Loop through s5255 and s52512 Results
        for index, scenario in s52512_csv_data.iloc[::-1].iterrows():

            # Get file name from scenario
            filename = scenario["File_Name"]

            if VERBOSE: print(f"\nFinding results for: {filename}")

            for oos_dir, is_dir in RESULT_DIRS:

                    plot_dir = os.path.join(output_dir, is_dir, "plots")
                    if not os.path.exists(plot_dir):
                        os.makedirs(plot_dir)

                    # Directories containing pickle files
                    oos_pkl_dir = os.path.join(output_dir, oos_dir, "pickle")
                    is_pkl_dir = os.path.join(output_dir, is_dir, "pickle")

                    # Locate pickle file that matches filename
                    nem_is_pkl_file = self.find_pickle(is_pkl_dir, filename)
                    nem_oos_pkl_file = self.find_pickle(oos_pkl_dir, filename)

                    # Skip if no pickle file found
                    if nem_is_pkl_file == None or nem_oos_pkl_file == None:
                        continue

                    if VERBOSE: print("\n", OUTPUT_DIR, is_dir, "plots", nem_is_pkl_file[:-4] + "_test_plot")

                    # Load in service and out of service data
                    nem_is_df = pd.read_pickle(os.path.join(is_pkl_dir,nem_is_pkl_file))
                    nem_oos_df = pd.read_pickle(os.path.join(oos_pkl_dir,nem_oos_pkl_file))

                    test_id = nem_is_pkl_file.replace("PSSE_s5.2.5.12","")[:-4]
                    load_case = is_dir

                    s52512_file_name = f"PSSE_S52512_NEM{test_id}.pdf"
                    s5255_file_name = f"PSSE_S5255_NEM{test_id}.pdf"
                    grid_file_name = f"NEM_grid_plot{test_id}.pdf"

                    s52512_save_path = os.path.join(output_dir, is_dir, "plots", s52512_file_name)
                    s5255_save_path = os.path.join(output_dir, is_dir, "plots", s5255_file_name)
                    grid_save_path = os.path.join(output_dir, is_dir, "plots", grid_file_name)

                    # # Perform Analysis on Scenario
                    # nem_analysis(
                    #     scenario=scenario,
                    #     results_df=nem_is_df,
                    # )

                    # # Add Columns to Results if analysis produces new columns.
                    # for col in scenario.index:
                    #     if col not in s52512_analysis_df.columns:
                    #         s52512_analysis_df[col] = ""

                    # # Save Results.
                    # s52512_analysis_df.loc[index] = scenario

                    # Generate Plots

                    if VERBOSE: print("\nGenerating S52512 Plot...")
                    self.in_and_out_of_service_plot("S52512 NEM Plot", load_case, test_id, scenario,nem_is_df, nem_oos_df, s52512_save_path)

                    #if VERBOSE: print("\nGenerating S5255 Plot...")
                    #self.s5255_plot("S5255 NEM Plot", load_case, test_id, scenario,nem_is_df, s5255_save_path)
                    # self.grid_plot(scenario, nem_is_df, grid_save_path)


        # Save Analysis Results
        output_directory = os.path.dirname("./_results/")

        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        # s52512_analysis_df.to_csv("./_results/s52512_nem_analysis.csv", index=False)
        # s52513_analysis_df.to_csv("./_results/s52513_nem_analysis.csv", index=False)

            


    def in_and_out_of_service_plot(self, title, load_case, test_id, scenario, nem_is_df, nem_oos_df, save_path):

        num_pages = 2
        num_of_cols = 2
        num_of_rows = 4

        signal_layout = [
            [ BUS_VOLTAGE_SIGNALS_1, BUS_VOLTAGE_SIGNALS_2,],
            [ACTIVE_POWER_SIGNALS, SYNC_MACHINE_ANGLE_SIGNALS,],
        ]

        # Specify the output PDF file path 
        with PdfPages(save_path) as pdf:

            # Loop through each page
            for page_number, page_signals in enumerate(signal_layout, 1):

            
                fig = self.mpl_setup(40,30) 

                if page_number == 1:

                    gs_main = gridspec.GridSpec(8, 5, figure=fig)
                    
                    cols = [
                        gridspec.GridSpecFromSubplotSpec(num_of_rows, 1, subplot_spec=gs_main[1:-1,0:3]),
                        gridspec.GridSpecFromSubplotSpec(num_of_rows, 1, subplot_spec=gs_main[1:-1,3:5]),
                    ]

                else:

                    gs_main = gridspec.GridSpec(8, 2, figure=fig)
                    
                    cols = [
                        gridspec.GridSpecFromSubplotSpec(num_of_rows, 1, subplot_spec=gs_main[1:-1,0]),
                        gridspec.GridSpecFromSubplotSpec(num_of_rows, 1, subplot_spec=gs_main[1:-1,1]),
                    ]

                description_ax = fig.add_subplot(gs_main[0:1,0:5])
                test_name = str(save_path).split('\\')[-3] + ", " + str(save_path).split('\\')[-1]
                description_ax.set_title(f"\n{title}: {load_case}{test_id} Page {page_number} of {num_pages}", style = 'italic', fontsize='large', loc='left', y=1.1, pad=0, color=AXIS_COLOUR)
                description_ax.axis("off")


                # TABLE PLOTTING

                if scenario['Fault_Bus'] == scenario['From_Bus']:
                    fault_loc = "From Bus"
                else:
                    fault_loc = "To Bus"


                fault_type_table = {
                    1:'1PHG',
                    4:'2PHG',
                    7:'3PHG',
                }

                fault_type = fault_type_table[int(scenario['Fault_Type'])]


                table_data = [
                    ['Filename/Int. Ref.', f"{scenario['File_Name']}"],
                    ['Contingency Type', f"{scenario['Contingency_Type']}"],
                    ['Contingency Loc', f"{scenario['Contingency_Loc']}"],
                    ['Fault Type', f"{fault_type}"],
                    ['Fault End', f"{fault_loc}"],
                    ['Auto Reclose', f"{scenario['Reclose']}"],
                ]

                colWidths = [0.20, 0.60]

                description_ax.table(cellText=table_data, colWidths=colWidths,cellLoc='left', bbox=[0,0,0.5,1])

                def set_table_edge_color(ax, color):
                    table: Table
                    for child in ax.get_children():
                        if isinstance(child, Table):
                            table = child
                            break

                    cell: Cell
                    for cell in table.get_children():
                        if isinstance(cell, Cell):
                            # print(cell)
                            cell.set_edgecolor(AXIS_COLOUR)
                            cell.set_linewidth(0.75)
                            cell.get_text().set_color(AXIS_COLOUR)

                set_table_edge_color(description_ax, AXIS_COLOUR)


                # Windlab Logo
                #self.plot_logo(fig)


                # Loop through each column
                for col_id, col_signals in enumerate(page_signals):

                    # Loop through each row
                    for row_id, signal_name in enumerate(col_signals):
                        
                        signal_name = signal_name[0]

                        ax = fig.add_subplot(cols[col_id][row_id,-1])

                        # Grab and Isolate Signals
                        nem_is_signal_raw = nem_is_df[signal_name][::PSSE_DECIMATE]

                        if "POC" in signal_name:
                            nem_oos_t = []
                            nem_oos_v = []

                        else:
                            nem_oos_signal_raw = nem_oos_df[signal_name][::PSSE_DECIMATE]

                            nem_oos_t = nem_oos_signal_raw.index
                            nem_oos_v = nem_oos_signal_raw.values

                        nem_is_t = nem_is_signal_raw.index
                        nem_is_v = nem_is_signal_raw.values

                        # Wrap PSSE phase angle 
                        if "DEG" in signal_name:
                            nem_oos_v = self.wrap_angle(nem_oos_v)
                            nem_is_v = self.wrap_angle(nem_is_v)

                        # Plot Signals
                        ax.plot(nem_is_t, nem_is_v, lw=SIGNAL_LINE_WIDTH, c=COL_SIG_1, zorder=4, label="In Service")
                        ax.plot(nem_oos_t, nem_oos_v, lw=SIGNAL_LINE_WIDTH, c=COL_SIG_2, zorder=3, label="Out of Service")

                        xlims = [
                            int(str(scenario["x_axis_limits"]).split(",")[0]),
                            int(str(scenario["x_axis_limits"]).split(",")[1]),
                        ]

                        self.config_subplot(ax, xlims, legend_ncol=2, signal_name=signal_name)

                        # Plot Title
                        title = signal_name.replace("_"," ")
                        ax.set_title(title, style = 'italic', fontsize='x-small', loc='left', y=1.00, pad=5, color=AXIS_COLOUR)

                        # Active Power Recovery
                        if signal_name == "MT_FOX_275KV_V_PU":

                            ymin, ymax = ax.get_ylim()
                            if abs(ymax - ymin) < 0.02:
                                ymin = ymin - 0.015
                                ymax = ymax + +0.015

                            ax.set_ylim(ymin=ymin, ymax=ymax)
        
                pdf.savefig() 
                plt.close()    

            plt.close()

    def grid_plot(self, scenario, df, savefile_path, plot_start=0, plot_duration=None, test_specific_table=None,
                  decimate=DECIMATE, filename="Filename", load_case="none", test_id="none" ):

        """
        Function to generate the medium level of detail grid plot

        Signals to add:
            BESSn_I_pu
            WTn_FaultCode

        To DO:
            can we add our thresholds as dotted lines? turbine lvrt/hvrt in/out, bess avr/gov freeze, vmp freeze
            can we put the qrefs and prefs for the batteries at the front?


            filter on the poc Qref? it's not sampled or anything so it'll be really choppy

        """

        # line weight for reference signals.
        lw_norm = 1.2
        lw_ref = 0.707 * lw_norm

        # Colours for signals.
        cref, cpoc, cw1, cw2, cw3, cw4 = np.array([
            (221, 179, 16),
            (64, 83, 211),
            (181, 29, 20),
            (0, 190, 255),
            (0, 178, 93),
            (251, 73, 176),
        ]) / 255

        cvmp = cpoc
        cb1, cb2, cb3, cb4 = [cw1, cw2, cw3, cw4]
        cx1, cx2, cx3, cx4 = [cw1, cw2, cw3, cw4]

        # Page Size and Subplot Layout
        number_of_rows = 6

        fig = plt.figure(facecolor=PRIMARY_COLOUR)  # subplot(number_of_rows, 4)

        cm = 1 / 2.54
        fig.set_size_inches(42 * cm, 29.7 * cm)  # A4 Size Landscape.
        plt.subplots_adjust(left=0.04, right=0.975, bottom=0.05, top=0.96, wspace=0.12, hspace=0.4)

        gs_main = gridspec.GridSpec(1, 4, figure=fig)

        gs_col_1 = gridspec.GridSpecFromSubplotSpec(number_of_rows, 1, subplot_spec=gs_main[0])
        gs_col_2 = gridspec.GridSpecFromSubplotSpec(number_of_rows, 1, subplot_spec=gs_main[1])
        gs_col_3 = gridspec.GridSpecFromSubplotSpec(number_of_rows, 1, subplot_spec=gs_main[2])
        gs_col_4 = gridspec.GridSpecFromSubplotSpec(number_of_rows, 1, subplot_spec=gs_main[3])
        gs_other = gridspec.GridSpecFromSubplotSpec(6, 1, subplot_spec=gs_col_4[:number_of_rows - 1])
        # plot_duration = df["POC_V_PU"].index[-1]

        ax_poc = []
        for i in range(number_of_rows):
            ax_poc.append(fig.add_subplot(gs_col_1[i, -1]))

        ax_wtg = []
        for i in range(number_of_rows):
            ax_wtg.append(fig.add_subplot(gs_col_2[i, -1]))

        ax_bess = []
        for i in range(number_of_rows):
            ax_bess.append(fig.add_subplot(gs_col_3[i, -1]))

        ax_table_1 = fig.add_subplot(gs_other[0, -1])
        ax_table_2 = fig.add_subplot(gs_other[1, -1])
        # ax_bess_rms_current = fig.add_subplot(gs_other[2, -1])
        ax_oltc = fig.add_subplot(gs_other[2, -1])
        ax_tap_pos = fig.add_subplot(gs_other[3, -1])
        ax_frt_1 = fig.add_subplot(gs_other[4, -1])
        ax_frt_2 = fig.add_subplot(gs_other[5, -1])

        ax_frequency = fig.add_subplot(gs_col_4[-1, -1])

        # Disable Unused Plots
        ax_table_1.axis('off')
        ax_table_2.axis('off')
        ax_frt_1.axis('off')
        ax_frt_2.axis('off')

        # Column Titles
        ax_poc[0].set_title('POC', style='italic', fontsize='medium', loc='center', y=1.1, pad=5, color=AXIS_COLOUR)
        ax_wtg[0].set_title('WTGs', style='italic', fontsize='medium', loc='center', y=1.1, pad=5, color=AXIS_COLOUR)
        ax_bess[0].set_title('BESS', style='italic', fontsize='medium', loc='center', y=1.1, pad=5, color=AXIS_COLOUR)
        ax_table_1.set_title('Test Parameters', style='italic', fontsize='medium', loc='center', y=1.04, pad=5,
                             color=AXIS_COLOUR)
        # ax_table_2.set_title('Test Specific Parameters', style = 'italic', fontsize='medium', loc='center', y=1.04, pad=5, color=AXIS_COLOUR)
        # ax_bess_rms_current.set_title('Other Plots', style = 'italic', fontsize='medium', loc='center', y=1.2, pad=5, color=AXIS_COLOUR)

        # Get Plot x-axis Parameters
        if plot_duration:
            plot_end = plot_start + plot_duration
        else:
            plot_end = df.index[-1]
            plot_duration = plot_end - plot_start

        if df.index[-1] < plot_end:
            plot_end = df.index[-1]
            plot_duration = plot_end - plot_start

        if scenario['Fault_Bus'] == scenario['From_Bus']:
            fault_loc = "From Bus"
        else:
            fault_loc = "To Bus"

        # TABLE PLOTTING
        table_data = [
            ['Filename/Int. Ref.', f"{scenario['File_Name']}"],
            ['Contingency Type', f"{scenario['Contingency_Type']}"],
            ['Contingency Loc', f"{scenario['Contingency_Loc']}"],
            ['Fault Type', f"{scenario['Fault_Type']}"],
            ['Fault End', f"{fault_loc}"],
            ['Auto Reclose', f"{scenario['Reclose']}"],

        ]

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

        for signal, label in optional_table_information:
            if signal in df.columns:
                table_data.append([label, f"{df[signal].values[0]}"])

        table_1_colour = np.empty_like(table_data, dtype='object')
        for i, _ in enumerate(table_1_colour):
            table_1_colour[i] = [SECONDARY_COLOUR, SECONDARY_COLOUR]

        colWidths = [0.35, 0.65]
        ax_table_1.table(cellText=table_data, colWidths=colWidths, cellLoc='left', loc='upper right',
                         cellColours=table_1_colour)


        def set_table_edge_color(ax, color):
            table: Table
            for child in ax.get_children():
                if isinstance(child, Table):
                    table = child
                    break

            cell: Cell
            for cell in table.get_children():
                if isinstance(cell, Cell):
                    # print(cell)
                    cell.set_edgecolor(AXIS_COLOUR)
                    cell.set_linewidth(0.75)
                    cell.get_text().set_color(AXIS_COLOUR)

        set_table_edge_color(ax_table_1, AXIS_COLOUR)
        if test_specific_table:
            set_table_edge_color(ax_table_2, AXIS_COLOUR)

        # PPOC PLOTTING
        pref_mw = df["PREF_POC_MW"][plot_start:plot_end][::decimate]
        ppoc_mw = df['POC_P_MW'][plot_start:plot_end][::decimate]
        self.signal_plot(
            ax=ax_poc[0],
            title='POC: P [MW]',
            traces=[('Pref', lw_ref, cref, 1, pref_mw),
                    ('Ppoc', lw_norm, cpoc, 2, ppoc_mw), ],
            start=plot_start,
            duration=plot_duration,
            fixed_y_lims=self.filtered_y_limits(ppoc_mw, pref_mw, min_y_range=40),
        )

        # QPOC PLOTTING
        qpoc_mvar = df['POC_Q_MVAR'][plot_start:plot_end][::decimate]
        # if scenario['Control_Mode'] == QCtrlMode.VAR:
        #     qref_mvar = df['QREF_MVAR'][plot_start:][::decimate]
        # else:
        #     qref_mvar = df['QREF_DROOP_MVAR'][plot_start:][::decimate]
        qref_mvar = df['QREF_DROOP_MVAR'][plot_start:][::decimate]
        self.signal_plot(
            ax=ax_poc[1],
            title='POC: Q [MVAr]',
            traces=[('Qref', lw_ref, cref, 1, qref_mvar),
                    ('Qpoc', lw_norm, cpoc, 2, qpoc_mvar), ],
            start=plot_start,
            duration=plot_duration,
            fixed_y_lims=self.filtered_y_limits(qpoc_mvar, qref_mvar, min_y_range=20),
        )

        # VPOC PLOTTING
        vpoc_pu = df['POC_V_PU'][plot_start:plot_end][::decimate]
        vref_pu = df['VREF_PU'][plot_start:plot_end][::decimate]

        # ax_poc[2].hlines(UBWFNetwork.tuning['bess_lvrt_th_in'], 0, plot_duration, linestyle='-', lw=0.6, color='0.5',
        #                  zorder=1, label="BESS-Frz-in")
        # ax_poc[2].hlines(UBWFNetwork.tuning['bess_hvrt_th_in'], 0, plot_duration, linestyle='-', lw=0.6, color='0.5',
        #                  zorder=1)
        #
        # ax_poc[2].hlines(UBWFNetwork.tuning['bess_lvrt_th_out'], 0, plot_duration, linestyle='-.', lw=0.6, color='0.5',
        #                  zorder=1, label="BESS-Frz-out")
        # ax_poc[2].hlines(UBWFNetwork.tuning['bess_hvrt_th_out'], 0, plot_duration, linestyle='-.', lw=0.6, color='0.5',
        #                  zorder=1)

        self.signal_plot(
            ax=ax_poc[2],
            title='POC: V [.pu]',
            traces=[('Vref', lw_ref, cref, 1, vref_pu),
                    ('Vpoc', lw_norm, cpoc, 2, vpoc_pu), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # POC ID PLOTTING
        id_poc_pu = df['POC_ID_PU'][plot_start:plot_end][::decimate]
        self.signal_plot(
            ax=ax_poc[number_of_rows - 3],
            title='POC: Id [.pu]',
            traces=[('POC', lw_norm, cpoc, 1, id_poc_pu)],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # POC: Iq
        iq_poc_pu = df['POC_IQ_PU'][plot_start:plot_end][::decimate]
        self.signal_plot(
            ax=ax_poc[number_of_rows - 2],
            title='POC: Iq [.pu]',
            traces=[('POC', lw_norm, cpoc, 1, iq_poc_pu)],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # WTG P PLOTTING
        self.signal_plot(
            ax=ax_wtg[0],
            title='WTG: P [MW]',
            traces=[('Pref', lw_ref, cref, 5, df['PREF_WT_MW'][plot_start:plot_end][::decimate]),
                    ('W1', lw_norm, cw1, 4, df['WT1_P_MW'][plot_start:plot_end][::decimate]),
                    ('W2', lw_norm, cw2, 3, df['WT2_P_MW'][plot_start:plot_end][::decimate]),
                    ('W3', lw_norm, cw3, 2, df['WT3_P_MW'][plot_start:plot_end][::decimate]),
                    ('W4', lw_norm, cw4, 1, df['WT4_P_MW'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # WTG Q PLOTTING
        self.signal_plot(
            ax=ax_wtg[1],
            title='WTG: Q [MVAr]',
            traces=[('Qref', lw_ref, cref, 5, df['QREF_WT_MVAR'][plot_start:plot_end][::decimate]),
                    ('W1', lw_norm, cw1, 4, df['WT1_Q_MVAR'][plot_start:plot_end][::decimate]),
                    ('W2', lw_norm, cw2, 3, df['WT2_Q_MVAR'][plot_start:plot_end][::decimate]),
                    ('W3', lw_norm, cw3, 2, df['WT3_Q_MVAR'][plot_start:plot_end][::decimate]),
                    ('W4', lw_norm, cw4, 1, df['WT4_Q_MVAR'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # WTG V PLOTTING
        # ax_wtg[2].hlines(UBWFNetwork.tuning['lvrt_th_in'], 0, plot_duration, linestyle='-', lw=0.5, color='0.5',
        #                  zorder=1,
        #                  label="FRT-in")
        # ax_wtg[2].hlines(UBWFNetwork.tuning['hvrt_th_in'], 0, plot_duration, linestyle='-', lw=0.5, color='0.5',
        #                  zorder=1)
        # ax_wtg[2].hlines(UBWFNetwork.tuning['lvrt_th_out'], 0, plot_duration, linestyle='-.', lw=0.5, color='0.5',
        #                  zorder=1)
        # ax_wtg[2].hlines(UBWFNetwork.tuning['hvrt_th_out'], 0, plot_duration, linestyle='-.', lw=0.5, color='0.5',
        #                  zorder=1,
        #                  label="FRT-out")

        self.signal_plot(
            ax=ax_wtg[2],
            title='WTG: V [.pu]',
            traces=[('W1', lw_norm, cw1, 4, df['WT1_V_PU'][plot_start:plot_end][::decimate]),
                    ('W2', lw_norm, cw2, 3, df['WT2_V_PU'][plot_start:plot_end][::decimate]),
                    ('W3', lw_norm, cw3, 2, df['WT3_V_PU'][plot_start:plot_end][::decimate]),
                    ('W4', lw_norm, cw4, 1, df['WT4_V_PU'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # WTG Id PLOTTING
        self.signal_plot(
            ax=ax_wtg[number_of_rows - 3],
            title='WTG: Id [.pu]',
            traces=[('W1', lw_norm, cw1, 4, df['WT1_ID_PU'][plot_start:plot_end][::decimate]),
                    ('W2', lw_norm, cw2, 3, df['WT2_ID_PU'][plot_start:plot_end][::decimate]),
                    ('W3', lw_norm, cw3, 2, df['WT3_ID_PU'][plot_start:plot_end][::decimate]),
                    ('W4', lw_norm, cw4, 1, df['WT4_ID_PU'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # WTG Iq PLOTTING
        self.signal_plot(
            ax=ax_wtg[number_of_rows - 2],
            title='WTG: Iq [.pu]',
            traces=[('W1', lw_norm, cw1, 4, df['WT1_IQ_PU'][plot_start:plot_end][::decimate]),
                    ('W2', lw_norm, cw2, 3, df['WT2_IQ_PU'][plot_start:plot_end][::decimate]),
                    ('W3', lw_norm, cw3, 2, df['WT3_IQ_PU'][plot_start:plot_end][::decimate]),
                    ('W4', lw_norm, cw4, 1, df['WT4_IQ_PU'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # BESS P PLOTTING
        self.signal_plot(
            ax=ax_bess[0],
            title='BESS: P [MW]',
            traces=[('Pref', lw_ref, cref, 5, df["PREF_BESS_EA_MW"][plot_start:plot_end][::decimate]),
                    ('B1', lw_norm, cb1, 4, df['BESS1_P_MW'][plot_start:plot_end][::decimate]),
                    ('B2', lw_norm, cb2, 3, df['BESS2_P_MW'][plot_start:plot_end][::decimate]),
                    ('B3', lw_norm, cb3, 2, df['BESS3_P_MW'][plot_start:plot_end][::decimate]),
                    ('B4', lw_norm, cb4, 1, df['BESS4_P_MW'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.5,
        )

        # BESS Q PLOTTING

        self.signal_plot(
            ax=ax_bess[1],
            title='BESS: Q [MVAr]',
            traces=[
                # ('Qref', lw_ref, cref, 5, df["QREF_BESS_EA_MVAR"][plot_start:plot_end][::decimate]),
                ('B1', lw_norm, cb1, 4, df['BESS1_Q_MVAR'][plot_start:plot_end][::decimate]),
                ('B2', lw_norm, cb2, 3, df['BESS2_Q_MVAR'][plot_start:plot_end][::decimate]),
                ('B3', lw_norm, cb3, 2, df['BESS3_Q_MVAR'][plot_start:plot_end][::decimate]),
                ('B4', lw_norm, cb4, 1, df['BESS4_Q_MVAR'][plot_start:plot_end][::decimate]),
            ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=10,
        )

        # BESS V PLOTTING

        self.signal_plot(
            ax=ax_bess[2],
            title='BESS: V [.pu]',
            traces=[('B1', lw_norm, cb1, 4, df['BESS1_V_PU'][plot_start:plot_end][::decimate]),
                    ('B2', lw_norm, cb2, 3, df['BESS2_V_PU'][plot_start:plot_end][::decimate]),
                    ('B3', lw_norm, cb3, 2, df['BESS3_V_PU'][plot_start:plot_end][::decimate]),
                    ('B4', lw_norm, cb4, 1, df['BESS4_V_PU'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # BESS Id PLOTTING
        self.signal_plot(
            ax=ax_bess[number_of_rows - 3],
            title='BESS: Id [.pu]',
            traces=[('B1', lw_norm, cb1, 4, df['BESS1_ID_PU'][plot_start:plot_end][::decimate]),
                    ('B2', lw_norm, cb2, 3, df['BESS2_ID_PU'][plot_start:plot_end][::decimate]),
                    ('B3', lw_norm, cb3, 2, df['BESS3_ID_PU'][plot_start:plot_end][::decimate]),
                    ('B4', lw_norm, cb4, 1, df['BESS4_ID_PU'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # BESS Iq PLOTTING
        self.signal_plot(
            ax=ax_bess[number_of_rows - 2],
            title='BESS: Iq [.pu]',
            traces=[('B1', lw_norm, cb1, 4, df['BESS1_IQ_PU'][plot_start:plot_end][::decimate]),
                    ('B2', lw_norm, cb2, 3, df['BESS2_IQ_PU'][plot_start:plot_end][::decimate]),
                    ('B3', lw_norm, cb3, 2, df['BESS3_IQ_PU'][plot_start:plot_end][::decimate]),
                    ('B4', lw_norm, cb4, 1, df['BESS4_IQ_PU'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # ANGLE PLOTS
        self.signal_plot(
            ax=ax_poc[-1],
            title='POC: Angle [deg]',
            traces=[
                ("POC", lw_norm, cpoc, 5, self.wrap_angle(df['POC_ANG_DEG'][plot_start:plot_end][::decimate])),
                # ("SLACK", lw_norm, (0,1,0), 5, wrap_angle(df['SLACK_ANG_DEG'][plot_start:plot_end][::decimate])),
            ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.1,
            time_axis_on=True,
        )

        self.signal_plot(
            ax=ax_wtg[-1],
            title='WTG: Angle [deg]',
            traces=[('W1', lw_norm, cw1, 4, self.wrap_angle(df['WT1_ANG_DEG'][plot_start:plot_end][::decimate])),
                    ('W2', lw_norm, cw2, 3, self.wrap_angle(df['WT2_ANG_DEG'][plot_start:plot_end][::decimate])),
                    ('W3', lw_norm, cw3, 2, self.wrap_angle(df['WT3_ANG_DEG'][plot_start:plot_end][::decimate])),
                    ('W4', lw_norm, cw4, 1, self.wrap_angle(df['WT4_ANG_DEG'][plot_start:plot_end][::decimate])), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.1,
            time_axis_on=True,
        )
        self.signal_plot(
            ax=ax_bess[-1],
            title='BESS: Angle [deg]',
            traces=[('B1', lw_norm, cb1, 4, self.wrap_angle(df['BESS1_ANG_DEG'][plot_start:plot_end][::decimate])),
                    ('B2', lw_norm, cb2, 3, self.wrap_angle(df['BESS2_ANG_DEG'][plot_start:plot_end][::decimate])),
                    ('B3', lw_norm, cb3, 2, self.wrap_angle(df['BESS3_ANG_DEG'][plot_start:plot_end][::decimate])),
                    ('B4', lw_norm, cb4, 1, self.wrap_angle(df['BESS4_ANG_DEG'][plot_start:plot_end][::decimate])), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.1,
            time_axis_on=True,
        )

        self.signal_plot(
            ax=ax_frequency,
            title='POC: Frequency [Hz]',
            traces=[("POC", lw_norm, cpoc, 5, df['FPOC_HZ'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.5,
            time_axis_on=True,
        )

        # FRT & Trip Sub Plot

        # Merge two subplots
        gs = ax_frt_1.get_gridspec()
        ax_frt_1.remove()
        ax_frt_2.remove()
        ax_frt = fig.add_subplot(gs[4:6, -1])

        # FRT Plotting Parameters
        frt_decimate = 1  # High detail required to capture frt flags
        major_spacing = 11
        group_offset = 1.5  # Space between BESS and WTG plots
        frz_flag_color = (0.6, 0.6, 1)
        hvrt_flag_color = (1, 0.2, 0.2)
        lvrt_flag_color = (0.2, 0.2, 1)
        trip_flag_color = (0, 0, 0)
        center_line_color = (0.8, 0.8, 0.8)

        # Signal Processing

        # vmp_frt_mode = df['EMP_V_FRZ'][plot_start:plot_end][::frt_decimate]
        vmp_frt_mode = df['VMP_FRT'][plot_start:plot_end][::frt_decimate]

        w1_frt_mode = df['WT1_FRT'][plot_start:plot_end][::frt_decimate]
        w2_frt_mode = df['WT2_FRT'][plot_start:plot_end][::frt_decimate]
        w3_frt_mode = df['WT3_FRT'][plot_start:plot_end][::frt_decimate]
        w4_frt_mode = df['WT4_FRT'][plot_start:plot_end][::frt_decimate]

        w1_trip_signal = df['WT1_TRIP'][plot_start:plot_end][::frt_decimate]
        w2_trip_signal = df['WT2_TRIP'][plot_start:plot_end][::frt_decimate]
        w3_trip_signal = df['WT3_TRIP'][plot_start:plot_end][::frt_decimate]
        w4_trip_signal = df['WT4_TRIP'][plot_start:plot_end][::frt_decimate]

        b1_trip_signal = df['BESS1_TRIP'][plot_start:plot_end][::frt_decimate]
        b2_trip_signal = df['BESS2_TRIP'][plot_start:plot_end][::frt_decimate]
        b3_trip_signal = df['BESS3_TRIP'][plot_start:plot_end][::frt_decimate]
        b4_trip_signal = df['BESS4_TRIP'][plot_start:plot_end][::frt_decimate]

        _, vmp_frz = self.read_flag_signal(vmp_frt_mode)

        w1_lvrt, w1_hvrt = self.read_psse_frt_flags(w1_frt_mode)
        w2_lvrt, w2_hvrt = self.read_psse_frt_flags(w2_frt_mode)
        w3_lvrt, w3_hvrt = self.read_psse_frt_flags(w3_frt_mode)
        w4_lvrt, w4_hvrt = self.read_psse_frt_flags(w4_frt_mode)

        _, w1_trip = self.read_flag_signal(w1_trip_signal)
        _, w2_trip = self.read_flag_signal(w2_trip_signal)
        _, w3_trip = self.read_flag_signal(w3_trip_signal)
        _, w4_trip = self.read_flag_signal(w4_trip_signal)

        _, b1_trip = self.read_flag_signal(b1_trip_signal)
        _, b2_trip = self.read_flag_signal(b2_trip_signal)
        _, b3_trip = self.read_flag_signal(b3_trip_signal)
        _, b4_trip = self.read_flag_signal(b4_trip_signal)

        _, b1_frz = self.read_flag_signal(df['BESS1_V_FRZ'][plot_start:plot_end][::frt_decimate])
        _, b2_frz = self.read_flag_signal(df['BESS2_V_FRZ'][plot_start:plot_end][::frt_decimate])
        _, b3_frz = self.read_flag_signal(df['BESS3_V_FRZ'][plot_start:plot_end][::frt_decimate])
        _, b4_frz = self.read_flag_signal(df['BESS4_V_FRZ'][plot_start:plot_end][::frt_decimate])

        traces = [
            ('FRZ', frz_flag_color, [b4_frz, None, b3_frz, None, b2_frz, None, b1_frz, None, vmp_frz]),
            ('HVRT', hvrt_flag_color, [None, w4_hvrt, None, w3_hvrt, None, w2_hvrt, None, w1_hvrt, None]),
            ('LVRT', lvrt_flag_color, [None, w4_lvrt, None, w3_lvrt, None, w2_lvrt, None, w1_lvrt, None]),
            ('Trip', trip_flag_color,
             [None, w4_trip, b4_trip, w3_trip, b3_trip, w2_trip, b2_trip, w1_trip, b1_trip, None]),
        ]

        # Centerline and Label Positions
        label_positions = []
        y_positions = []

        for y in range(1, 5):
            y_positions.append(y * major_spacing - group_offset)
            # y_positions.append(y * major_spacing )
            y_positions.append(y * major_spacing + group_offset)
            label_positions.append(y * major_spacing)

        y_positions.append(5 * major_spacing - group_offset)
        label_positions.append(5 * major_spacing - group_offset)

        ax_frt.hlines(y_positions, 0, plot_duration, linestyle='--', lw=0.45, color=center_line_color, zorder=1)

        # Plot Traces
        for (label, flag_color, flags) in traces:

            ax_frt.plot([], [], c=flag_color, label=label)

            for index, flag in enumerate(flags):
                if flag is not None:
                    for flag_segment in flag:
                        self.plot_flag(ax_frt, plot_start, flag_segment, y_positions[index], flag_color)

        self.configure_subplot(ax_frt, 'VMP/WTG/BESS FRT & Trip Flags', plot_duration,
                          fixed_y_lims=(group_offset, 6 * major_spacing - 3 * group_offset))

        ax_frt.set_yticks(y_positions, ['B4', 'W4', 'B3', 'W3', 'B2', 'W2', 'B1', 'W1', 'VMP'])
        ax_frt.grid(axis='y')
        ax_frt.tick_params(axis='y', labelsize=6.0, left=False)

        traces = [('TX1', cx1, 4, df['TX1_TAP'][plot_start:plot_end][::decimate]),
                  ('TX2', cx2, 3, df['TX2_TAP'][plot_start:plot_end][::decimate]),
                  ('TX3', cx3, 2, df['TX3_TAP'][plot_start:plot_end][::decimate]),
                  ('TX4', cx4, 1, df['TX4_TAP'][plot_start:plot_end][::decimate]),
                  ]
        for (label, tcolor, zorder, signal) in traces:
            ax_tap_pos.plot(signal.index - plot_start, signal.values,
                            lw=lw_norm, c=tcolor, zorder=zorder, label=label)
        self.configure_subplot(ax_tap_pos, 'TX Taps Positions', plot_duration,
                          legendon=True, xlabelson=False, min_y_range=0.05)

        # OLTC Tapping

        # Merge two subplots
        gs_tx = ax_oltc.get_gridspec()
        # ax_tap_pos.remove()
        ax_oltc.remove()
        ax_tx = fig.add_subplot(gs_tx[2:3, -1])

        # FRT Plotting Parameters
        oltc_decimate = 5  # High detail required to capture frt flags
        major_spacing = 10
        group_offset = 3  # Space between BESS and WTG plots
        # tapping_flag_color = (0.9,0.2,0.2)
        tapping_flag_color = '0.65'  # light grey
        tap_pos_color = (0.9, 0.9, 0.9)
        center_line_color = (0.8, 0.8, 0.8)

        # Signal Processing

        _, tx1_settled = self.read_flag_signal(df['TX1_SETTLED'][plot_start:plot_end][::oltc_decimate])
        _, tx2_settled = self.read_flag_signal(df['TX2_SETTLED'][plot_start:plot_end][::oltc_decimate])
        _, tx3_settled = self.read_flag_signal(df['TX3_SETTLED'][plot_start:plot_end][::oltc_decimate])
        _, tx4_settled = self.read_flag_signal(df['TX4_SETTLED'][plot_start:plot_end][::oltc_decimate])

        traces = [
            ('TX Settled', tapping_flag_color, [tx1_settled, tx2_settled, tx3_settled, tx4_settled]),
        ]

        # Centerline and Label Positions
        label_positions = []
        y_positions = []

        for y in range(1, 5):
            # y_positions.append(y * major_spacing - group_offset)
            y_positions.append(y * major_spacing)
            # y_positions.append(y * major_spacing + group_offset)
            label_positions.append(y * major_spacing)

        # y_positions.append(5*major_spacing - group_offset)
        # label_positions.append(5*major_spacing - group_offset)

        ax_tx.hlines(y_positions, 0, plot_duration, linestyle='--', lw=0.6, color=center_line_color, zorder=1)

        # Plot Traces
        for (label, flag_color, flags) in traces:

            ax_tx.plot([], [], c=flag_color, label=label)

            for index, flag in enumerate(flags):
                if flag is not None:
                    for flag_segment in flag:
                        self.plot_flag(ax_tx, plot_start, flag_segment, y_positions[index], flag_color)

        self.configure_subplot(ax_tx, 'OLTC Tapping Flags', plot_duration,
                          fixed_y_lims=(group_offset, 5 * major_spacing - group_offset))

        ax_tx.set_yticks(y_positions, ['TX4', 'TX3', 'TX2', 'TX1'])
        ax_tx.grid(axis='y')
        ax_tx.tick_params(axis='y', labelsize=6.0, left=False)

        fig.align_labels()

        plt.savefig(savefile_path, bbox_inches='tight', dpi=150, format='pdf')
        # plt.savefig(savefile_path + '.pdf', bbox_inches='tight', dpi=150, format='pdf')
        # plt.savefig(savefile_path + '.png', bbox_inches='tight', dpi=150, format='png')
        plt.cla()  # Clear the current axes.
        plt.clf()  # Clear the current figure.
        plt.close()  # Closes all the figure windows.

        """
        Function to generate the medium level of detail grid plot

        Signals to add:
            BESSn_I_pu
            WTn_FaultCode

        To DO:
            can we add our thresholds as dotted lines? turbine lvrt/hvrt in/out, bess avr/gov freeze, vmp freeze
            can we put the qrefs and prefs for the batteries at the front?


            filter on the poc Qref? it's not sampled or anything so it'll be really choppy

        """

        # line weight for reference signals.
        lw_norm = 1.2
        lw_ref = 0.707 * lw_norm

        # Colours for signals.
        cref, cpoc, cw1, cw2, cw3, cw4 = np.array([
            (221, 179, 16),
            (64, 83, 211),
            (181, 29, 20),
            (0, 190, 255),
            (0, 178, 93),
            (251, 73, 176),
        ]) / 255

        cvmp = cpoc
        cb1, cb2, cb3, cb4 = [cw1, cw2, cw3, cw4]
        cx1, cx2, cx3, cx4 = [cw1, cw2, cw3, cw4]

        # Page Size and Subplot Layout
        number_of_rows = 6

        fig = plt.figure(facecolor=PRIMARY_COLOUR)  # subplot(number_of_rows, 4)

        cm = 1 / 2.54
        fig.set_size_inches(42 * cm, 29.7 * cm)  # A4 Size Landscape.
        plt.subplots_adjust(left=0.04, right=0.975, bottom=0.05, top=0.96, wspace=0.12, hspace=0.4)

        gs_main = gridspec.GridSpec(1, 4, figure=fig)

        gs_col_1 = gridspec.GridSpecFromSubplotSpec(number_of_rows, 1, subplot_spec=gs_main[0])
        gs_col_2 = gridspec.GridSpecFromSubplotSpec(number_of_rows, 1, subplot_spec=gs_main[1])
        gs_col_3 = gridspec.GridSpecFromSubplotSpec(number_of_rows, 1, subplot_spec=gs_main[2])
        gs_col_4 = gridspec.GridSpecFromSubplotSpec(number_of_rows, 1, subplot_spec=gs_main[3])
        gs_other = gridspec.GridSpecFromSubplotSpec(6, 1, subplot_spec=gs_col_4[:number_of_rows - 1])
        # plot_duration = df["POC_V_PU"].index[-1]

        ax_poc = []
        for i in range(number_of_rows):
            ax_poc.append(fig.add_subplot(gs_col_1[i, -1]))

        ax_wtg = []
        for i in range(number_of_rows):
            ax_wtg.append(fig.add_subplot(gs_col_2[i, -1]))

        ax_bess = []
        for i in range(number_of_rows):
            ax_bess.append(fig.add_subplot(gs_col_3[i, -1]))

        ax_table_1 = fig.add_subplot(gs_other[0, -1])
        ax_table_2 = fig.add_subplot(gs_other[1, -1])
        # ax_bess_rms_current = fig.add_subplot(gs_other[2, -1])
        ax_oltc = fig.add_subplot(gs_other[2, -1])
        ax_tap_pos = fig.add_subplot(gs_other[3, -1])
        ax_frt_1 = fig.add_subplot(gs_other[4, -1])
        ax_frt_2 = fig.add_subplot(gs_other[5, -1])

        ax_frequency = fig.add_subplot(gs_col_4[-1, -1])

        # Disable Unused Plots
        ax_table_1.axis('off')
        ax_table_2.axis('off')
        ax_frt_1.axis('off')
        ax_frt_2.axis('off')

        # Column Titles
        ax_poc[0].set_title('POC', style='italic', fontsize='medium', loc='center', y=1.1, pad=5, color=AXIS_COLOUR)
        ax_wtg[0].set_title('WTGs', style='italic', fontsize='medium', loc='center', y=1.1, pad=5, color=AXIS_COLOUR)
        ax_bess[0].set_title('BESS', style='italic', fontsize='medium', loc='center', y=1.1, pad=5, color=AXIS_COLOUR)
        ax_table_1.set_title('Test Parameters', style='italic', fontsize='medium', loc='center', y=1.04, pad=5,
                             color=AXIS_COLOUR)
        # ax_table_2.set_title('Test Specific Parameters', style = 'italic', fontsize='medium', loc='center', y=1.04, pad=5, color=AXIS_COLOUR)
        # ax_bess_rms_current.set_title('Other Plots', style = 'italic', fontsize='medium', loc='center', y=1.2, pad=5, color=AXIS_COLOUR)

        # Get Plot x-axis Parameters
        if plot_duration:
            plot_end = plot_start + plot_duration
        else:
            plot_end = df.index[-1]
            plot_duration = plot_end - plot_start

        if df.index[-1] < plot_end:
            plot_end = df.index[-1]
            plot_duration = plot_end - plot_start

        # TABLE PLOTTING
        table_data = [
            ['Filename/Int. Ref.', f"{scenario['File_Name']}"],
            # ['Fault Level [MVA]', f"{scenario['Init_Fault_MVA']}"],
            # ['X2R', f"{scenario['Init_Fault_X_on_R']}"],
            # ['Vpoc [.pu]', f"{scenario['Init_Vpoc_pu']:.4f}"],
            # ['Qpoc [.pu]', f"{scenario['Init_Qpoc_pu']}"],
            # ['Pwind [MW]', f"{scenario['Init_Pwind_MW']}"],
            # ['Pbess [MW]', f"{scenario['Init_Pbess_MW']}"],
        ]

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

        for signal, label in optional_table_information:
            if signal in df.columns:
                table_data.append([label, f"{df[signal].values[0]}"])

        table_1_colour = np.empty_like(table_data, dtype='object')
        for i, _ in enumerate(table_1_colour):
            table_1_colour[i] = [SECONDARY_COLOUR, SECONDARY_COLOUR]

        colWidths = [0.35, 0.65]
        ax_table_1.table(cellText=table_data, colWidths=colWidths, cellLoc='left', loc='upper right',
                         cellColours=table_1_colour)


        def set_table_edge_color(ax, color):
            table: Table
            for child in ax.get_children():
                if isinstance(child, Table):
                    table = child
                    break

            cell: Cell
            for cell in table.get_children():
                if isinstance(cell, Cell):
                    # print(cell)
                    cell.set_edgecolor(AXIS_COLOUR)
                    cell.set_linewidth(0.75)
                    cell.get_text().set_color(AXIS_COLOUR)

        set_table_edge_color(ax_table_1, AXIS_COLOUR)
        if test_specific_table:
            set_table_edge_color(ax_table_2, AXIS_COLOUR)

        # PPOC PLOTTING
        pref_mw = df["PREF_POC_MW"][plot_start:plot_end][::decimate]
        ppoc_mw = df['POC_P_MW'][plot_start:plot_end][::decimate]
        self.signal_plot(
            ax=ax_poc[0],
            title='POC: P [MW]',
            traces=[('Pref', lw_ref, cref, 1, pref_mw),
                    ('Ppoc', lw_norm, cpoc, 2, ppoc_mw), ],
            start=plot_start,
            duration=plot_duration,
            fixed_y_lims=self.filtered_y_limits(ppoc_mw, pref_mw, min_y_range=40),
        )

        # QPOC PLOTTING
        qpoc_mvar = df['POC_Q_MVAR'][plot_start:plot_end][::decimate]
        # if scenario['Control_Mode'] == QCtrlMode.VAR:
        #     qref_mvar = df['QREF_MVAR'][plot_start:][::decimate]
        # else:
        #     qref_mvar = df['QREF_DROOP_MVAR'][plot_start:][::decimate]
        qref_mvar = df['QREF_DROOP_MVAR'][plot_start:][::decimate]
        self.signal_plot(
            ax=ax_poc[1],
            title='POC: Q [MVAr]',
            traces=[('Qref', lw_ref, cref, 1, qref_mvar),
                    ('Qpoc', lw_norm, cpoc, 2, qpoc_mvar), ],
            start=plot_start,
            duration=plot_duration,
            fixed_y_lims=self.filtered_y_limits(qpoc_mvar, qref_mvar, min_y_range=20),
        )

        # VPOC PLOTTING
        vpoc_pu = df['POC_V_PU'][plot_start:plot_end][::decimate]
        vref_pu = df['VREF_PU'][plot_start:plot_end][::decimate]

        # ax_poc[2].hlines(UBWFNetwork.tuning['bess_lvrt_th_in'], 0, plot_duration, linestyle='-', lw=0.6, color='0.5',
        #                  zorder=1, label="BESS-Frz-in")
        # ax_poc[2].hlines(UBWFNetwork.tuning['bess_hvrt_th_in'], 0, plot_duration, linestyle='-', lw=0.6, color='0.5',
        #                  zorder=1)
        #
        # ax_poc[2].hlines(UBWFNetwork.tuning['bess_lvrt_th_out'], 0, plot_duration, linestyle='-.', lw=0.6, color='0.5',
        #                  zorder=1, label="BESS-Frz-out")
        # ax_poc[2].hlines(UBWFNetwork.tuning['bess_hvrt_th_out'], 0, plot_duration, linestyle='-.', lw=0.6, color='0.5',
        #                  zorder=1)

        self.signal_plot(
            ax=ax_poc[2],
            title='POC: V [.pu]',
            traces=[('Vref', lw_ref, cref, 1, vref_pu),
                    ('Vpoc', lw_norm, cpoc, 2, vpoc_pu), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # POC ID PLOTTING
        id_poc_pu = df['POC_ID_PU'][plot_start:plot_end][::decimate]
        self.signal_plot(
            ax=ax_poc[number_of_rows - 3],
            title='POC: Id [.pu]',
            traces=[('POC', lw_norm, cpoc, 1, id_poc_pu)],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # POC: Iq
        iq_poc_pu = df['POC_IQ_PU'][plot_start:plot_end][::decimate]
        self.signal_plot(
            ax=ax_poc[number_of_rows - 2],
            title='POC: Iq [.pu]',
            traces=[('POC', lw_norm, cpoc, 1, iq_poc_pu)],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # WTG P PLOTTING
        self.signal_plot(
            ax=ax_wtg[0],
            title='WTG: P [MW]',
            traces=[('Pref', lw_ref, cref, 5, df['PREF_WT_MW'][plot_start:plot_end][::decimate]),
                    ('W1', lw_norm, cw1, 4, df['WT1_P_MW'][plot_start:plot_end][::decimate]),
                    ('W2', lw_norm, cw2, 3, df['WT2_P_MW'][plot_start:plot_end][::decimate]),
                    ('W3', lw_norm, cw3, 2, df['WT3_P_MW'][plot_start:plot_end][::decimate]),
                    ('W4', lw_norm, cw4, 1, df['WT4_P_MW'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # WTG Q PLOTTING
        self.signal_plot(
            ax=ax_wtg[1],
            title='WTG: Q [MVAr]',
            traces=[('Qref', lw_ref, cref, 5, df['QREF_WT_MVAR'][plot_start:plot_end][::decimate]),
                    ('W1', lw_norm, cw1, 4, df['WT1_Q_MVAR'][plot_start:plot_end][::decimate]),
                    ('W2', lw_norm, cw2, 3, df['WT2_Q_MVAR'][plot_start:plot_end][::decimate]),
                    ('W3', lw_norm, cw3, 2, df['WT3_Q_MVAR'][plot_start:plot_end][::decimate]),
                    ('W4', lw_norm, cw4, 1, df['WT4_Q_MVAR'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # WTG V PLOTTING
        # ax_wtg[2].hlines(UBWFNetwork.tuning['lvrt_th_in'], 0, plot_duration, linestyle='-', lw=0.5, color='0.5',
        #                  zorder=1,
        #                  label="FRT-in")
        # ax_wtg[2].hlines(UBWFNetwork.tuning['hvrt_th_in'], 0, plot_duration, linestyle='-', lw=0.5, color='0.5',
        #                  zorder=1)
        # ax_wtg[2].hlines(UBWFNetwork.tuning['lvrt_th_out'], 0, plot_duration, linestyle='-.', lw=0.5, color='0.5',
        #                  zorder=1)
        # ax_wtg[2].hlines(UBWFNetwork.tuning['hvrt_th_out'], 0, plot_duration, linestyle='-.', lw=0.5, color='0.5',
        #                  zorder=1,
        #                  label="FRT-out")

        self.signal_plot(
            ax=ax_wtg[2],
            title='WTG: V [.pu]',
            traces=[('W1', lw_norm, cw1, 4, df['WT1_V_PU'][plot_start:plot_end][::decimate]),
                    ('W2', lw_norm, cw2, 3, df['WT2_V_PU'][plot_start:plot_end][::decimate]),
                    ('W3', lw_norm, cw3, 2, df['WT3_V_PU'][plot_start:plot_end][::decimate]),
                    ('W4', lw_norm, cw4, 1, df['WT4_V_PU'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # WTG Id PLOTTING
        self.signal_plot(
            ax=ax_wtg[number_of_rows - 3],
            title='WTG: Id [.pu]',
            traces=[('W1', lw_norm, cw1, 4, df['WT1_ID_PU'][plot_start:plot_end][::decimate]),
                    ('W2', lw_norm, cw2, 3, df['WT2_ID_PU'][plot_start:plot_end][::decimate]),
                    ('W3', lw_norm, cw3, 2, df['WT3_ID_PU'][plot_start:plot_end][::decimate]),
                    ('W4', lw_norm, cw4, 1, df['WT4_ID_PU'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # WTG Iq PLOTTING
        self.signal_plot(
            ax=ax_wtg[number_of_rows - 2],
            title='WTG: Iq [.pu]',
            traces=[('W1', lw_norm, cw1, 4, df['WT1_IQ_PU'][plot_start:plot_end][::decimate]),
                    ('W2', lw_norm, cw2, 3, df['WT2_IQ_PU'][plot_start:plot_end][::decimate]),
                    ('W3', lw_norm, cw3, 2, df['WT3_IQ_PU'][plot_start:plot_end][::decimate]),
                    ('W4', lw_norm, cw4, 1, df['WT4_IQ_PU'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # BESS P PLOTTING
        self.signal_plot(
            ax=ax_bess[0],
            title='BESS: P [MW]',
            traces=[('Pref', lw_ref, cref, 5, df["PREF_BESS_EA_MW"][plot_start:plot_end][::decimate]),
                    ('B1', lw_norm, cb1, 4, df['BESS1_P_MW'][plot_start:plot_end][::decimate]),
                    ('B2', lw_norm, cb2, 3, df['BESS2_P_MW'][plot_start:plot_end][::decimate]),
                    ('B3', lw_norm, cb3, 2, df['BESS3_P_MW'][plot_start:plot_end][::decimate]),
                    ('B4', lw_norm, cb4, 1, df['BESS4_P_MW'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.5,
        )

        # BESS Q PLOTTING

        self.signal_plot(
            ax=ax_bess[1],
            title='BESS: Q [MVAr]',
            traces=[
                # ('Qref', lw_ref, cref, 5, df["QREF_BESS_EA_MVAR"][plot_start:plot_end][::decimate]),
                ('B1', lw_norm, cb1, 4, df['BESS1_Q_MVAR'][plot_start:plot_end][::decimate]),
                ('B2', lw_norm, cb2, 3, df['BESS2_Q_MVAR'][plot_start:plot_end][::decimate]),
                ('B3', lw_norm, cb3, 2, df['BESS3_Q_MVAR'][plot_start:plot_end][::decimate]),
                ('B4', lw_norm, cb4, 1, df['BESS4_Q_MVAR'][plot_start:plot_end][::decimate]),
            ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=10,
        )

        # BESS V PLOTTING

        self.signal_plot(
            ax=ax_bess[2],
            title='BESS: V [.pu]',
            traces=[('B1', lw_norm, cb1, 4, df['BESS1_V_PU'][plot_start:plot_end][::decimate]),
                    ('B2', lw_norm, cb2, 3, df['BESS2_V_PU'][plot_start:plot_end][::decimate]),
                    ('B3', lw_norm, cb3, 2, df['BESS3_V_PU'][plot_start:plot_end][::decimate]),
                    ('B4', lw_norm, cb4, 1, df['BESS4_V_PU'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # BESS Id PLOTTING
        self.signal_plot(
            ax=ax_bess[number_of_rows - 3],
            title='BESS: Id [.pu]',
            traces=[('B1', lw_norm, cb1, 4, df['BESS1_ID_PU'][plot_start:plot_end][::decimate]),
                    ('B2', lw_norm, cb2, 3, df['BESS2_ID_PU'][plot_start:plot_end][::decimate]),
                    ('B3', lw_norm, cb3, 2, df['BESS3_ID_PU'][plot_start:plot_end][::decimate]),
                    ('B4', lw_norm, cb4, 1, df['BESS4_ID_PU'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # BESS Iq PLOTTING
        self.signal_plot(
            ax=ax_bess[number_of_rows - 2],
            title='BESS: Iq [.pu]',
            traces=[('B1', lw_norm, cb1, 4, df['BESS1_IQ_PU'][plot_start:plot_end][::decimate]),
                    ('B2', lw_norm, cb2, 3, df['BESS2_IQ_PU'][plot_start:plot_end][::decimate]),
                    ('B3', lw_norm, cb3, 2, df['BESS3_IQ_PU'][plot_start:plot_end][::decimate]),
                    ('B4', lw_norm, cb4, 1, df['BESS4_IQ_PU'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.05,
        )

        # ANGLE PLOTS
        self.signal_plot(
            ax=ax_poc[-1],
            title='POC: Angle [deg]',
            traces=[
                ("POC", lw_norm, cpoc, 5, self.wrap_angle(df['POC_ANG_DEG'][plot_start:plot_end][::decimate])),
                # ("SLACK", lw_norm, (0,1,0), 5, wrap_angle(df['SLACK_ANG_DEG'][plot_start:plot_end][::decimate])),
            ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.1,
            time_axis_on=True,
        )

        self.signal_plot(
            ax=ax_wtg[-1],
            title='WTG: Angle [deg]',
            traces=[('W1', lw_norm, cw1, 4, self.wrap_angle(df['WT1_ANG_DEG'][plot_start:plot_end][::decimate])),
                    ('W2', lw_norm, cw2, 3, self.wrap_angle(df['WT2_ANG_DEG'][plot_start:plot_end][::decimate])),
                    ('W3', lw_norm, cw3, 2, self.wrap_angle(df['WT3_ANG_DEG'][plot_start:plot_end][::decimate])),
                    ('W4', lw_norm, cw4, 1, self.wrap_angle(df['WT4_ANG_DEG'][plot_start:plot_end][::decimate])), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.1,
            time_axis_on=True,
        )
        self.signal_plot(
            ax=ax_bess[-1],
            title='BESS: Angle [deg]',
            traces=[('B1', lw_norm, cb1, 4, self.wrap_angle(df['BESS1_ANG_DEG'][plot_start:plot_end][::decimate])),
                    ('B2', lw_norm, cb2, 3, self.wrap_angle(df['BESS2_ANG_DEG'][plot_start:plot_end][::decimate])),
                    ('B3', lw_norm, cb3, 2, self.wrap_angle(df['BESS3_ANG_DEG'][plot_start:plot_end][::decimate])),
                    ('B4', lw_norm, cb4, 1, self.wrap_angle(df['BESS4_ANG_DEG'][plot_start:plot_end][::decimate])), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.1,
            time_axis_on=True,
        )

        self.signal_plot(
            ax=ax_frequency,
            title='POC: Frequency [Hz]',
            traces=[("POC", lw_norm, cpoc, 5, df['FPOC_HZ'][plot_start:plot_end][::decimate]), ],
            start=plot_start,
            duration=plot_duration,
            min_y_range=0.5,
            time_axis_on=True,
        )

        # FRT & Trip Sub Plot

        # Merge two subplots
        gs = ax_frt_1.get_gridspec()
        ax_frt_1.remove()
        ax_frt_2.remove()
        ax_frt = fig.add_subplot(gs[4:6, -1])

        # FRT Plotting Parameters
        frt_decimate = 1  # High detail required to capture frt flags
        major_spacing = 11
        group_offset = 1.5  # Space between BESS and WTG plots
        frz_flag_color = (0.6, 0.6, 1)
        hvrt_flag_color = (1, 0.2, 0.2)
        lvrt_flag_color = (0.2, 0.2, 1)
        trip_flag_color = (0, 0, 0)
        center_line_color = (0.8, 0.8, 0.8)

        # Signal Processing

        # vmp_frt_mode = df['EMP_V_FRZ'][plot_start:plot_end][::frt_decimate]
        vmp_frt_mode = df['VMP_FRT'][plot_start:plot_end][::frt_decimate]

        w1_frt_mode = df['WT1_FRT'][plot_start:plot_end][::frt_decimate]
        w2_frt_mode = df['WT2_FRT'][plot_start:plot_end][::frt_decimate]
        w3_frt_mode = df['WT3_FRT'][plot_start:plot_end][::frt_decimate]
        w4_frt_mode = df['WT4_FRT'][plot_start:plot_end][::frt_decimate]

        w1_trip_signal = df['WT1_TRIP'][plot_start:plot_end][::frt_decimate]
        w2_trip_signal = df['WT2_TRIP'][plot_start:plot_end][::frt_decimate]
        w3_trip_signal = df['WT3_TRIP'][plot_start:plot_end][::frt_decimate]
        w4_trip_signal = df['WT4_TRIP'][plot_start:plot_end][::frt_decimate]

        b1_trip_signal = df['BESS1_TRIP'][plot_start:plot_end][::frt_decimate]
        b2_trip_signal = df['BESS2_TRIP'][plot_start:plot_end][::frt_decimate]
        b3_trip_signal = df['BESS3_TRIP'][plot_start:plot_end][::frt_decimate]
        b4_trip_signal = df['BESS4_TRIP'][plot_start:plot_end][::frt_decimate]

        _, vmp_frz = self.read_flag_signal(vmp_frt_mode)

        w1_lvrt, w1_hvrt = self.read_psse_frt_flags(w1_frt_mode)
        w2_lvrt, w2_hvrt = self.read_psse_frt_flags(w2_frt_mode)
        w3_lvrt, w3_hvrt = self.read_psse_frt_flags(w3_frt_mode)
        w4_lvrt, w4_hvrt = self.read_psse_frt_flags(w4_frt_mode)

        _, w1_trip = self.read_flag_signal(w1_trip_signal)
        _, w2_trip = self.read_flag_signal(w2_trip_signal)
        _, w3_trip = self.read_flag_signal(w3_trip_signal)
        _, w4_trip = self.read_flag_signal(w4_trip_signal)

        _, b1_trip = self.read_flag_signal(b1_trip_signal)
        _, b2_trip = self.read_flag_signal(b2_trip_signal)
        _, b3_trip = self.read_flag_signal(b3_trip_signal)
        _, b4_trip = self.read_flag_signal(b4_trip_signal)

        _, b1_frz = self.read_flag_signal(df['BESS1_V_FRZ'][plot_start:plot_end][::frt_decimate])
        _, b2_frz = self.read_flag_signal(df['BESS2_V_FRZ'][plot_start:plot_end][::frt_decimate])
        _, b3_frz = self.read_flag_signal(df['BESS3_V_FRZ'][plot_start:plot_end][::frt_decimate])
        _, b4_frz = self.read_flag_signal(df['BESS4_V_FRZ'][plot_start:plot_end][::frt_decimate])

        traces = [
            ('FRZ', frz_flag_color, [b4_frz, None, b3_frz, None, b2_frz, None, b1_frz, None, vmp_frz]),
            ('HVRT', hvrt_flag_color, [None, w4_hvrt, None, w3_hvrt, None, w2_hvrt, None, w1_hvrt, None]),
            ('LVRT', lvrt_flag_color, [None, w4_lvrt, None, w3_lvrt, None, w2_lvrt, None, w1_lvrt, None]),
            ('Trip', trip_flag_color,
             [None, w4_trip, b4_trip, w3_trip, b3_trip, w2_trip, b2_trip, w1_trip, b1_trip, None]),
        ]

        # Centerline and Label Positions
        label_positions = []
        y_positions = []

        for y in range(1, 5):
            y_positions.append(y * major_spacing - group_offset)
            # y_positions.append(y * major_spacing )
            y_positions.append(y * major_spacing + group_offset)
            label_positions.append(y * major_spacing)

        y_positions.append(5 * major_spacing - group_offset)
        label_positions.append(5 * major_spacing - group_offset)

        ax_frt.hlines(y_positions, 0, plot_duration, linestyle='--', lw=0.45, color=center_line_color, zorder=1)

        # Plot Traces
        for (label, flag_color, flags) in traces:

            ax_frt.plot([], [], c=flag_color, label=label)

            for index, flag in enumerate(flags):
                if flag is not None:
                    for flag_segment in flag:
                        self.plot_flag(ax_frt, plot_start, flag_segment, y_positions[index], flag_color)

        self.configure_subplot(ax_frt, 'VMP/WTG/BESS FRT & Trip Flags', plot_duration,
                          fixed_y_lims=(group_offset, 6 * major_spacing - 3 * group_offset))

        ax_frt.set_yticks(y_positions, ['B4', 'W4', 'B3', 'W3', 'B2', 'W2', 'B1', 'W1', 'VMP'])
        ax_frt.grid(axis='y')
        ax_frt.tick_params(axis='y', labelsize=6.0, left=False)

        traces = [('TX1', cx1, 4, df['TX1_TAP'][plot_start:plot_end][::decimate]),
                  ('TX2', cx2, 3, df['TX2_TAP'][plot_start:plot_end][::decimate]),
                  ('TX3', cx3, 2, df['TX3_TAP'][plot_start:plot_end][::decimate]),
                  ('TX4', cx4, 1, df['TX4_TAP'][plot_start:plot_end][::decimate]),
                  ]
        for (label, tcolor, zorder, signal) in traces:
            ax_tap_pos.plot(signal.index - plot_start, signal.values,
                            lw=lw_norm, c=tcolor, zorder=zorder, label=label)
        self.configure_subplot(ax_tap_pos, 'TX Taps Positions', plot_duration,
                          legendon=True, xlabelson=False, min_y_range=0.05)

        # OLTC Tapping

        # Merge two subplots
        gs_tx = ax_oltc.get_gridspec()
        # ax_tap_pos.remove()
        ax_oltc.remove()
        ax_tx = fig.add_subplot(gs_tx[2:3, -1])

        # FRT Plotting Parameters
        oltc_decimate = 5  # High detail required to capture frt flags
        major_spacing = 10
        group_offset = 3  # Space between BESS and WTG plots
        # tapping_flag_color = (0.9,0.2,0.2)
        tapping_flag_color = '0.65'  # light grey
        tap_pos_color = (0.9, 0.9, 0.9)
        center_line_color = (0.8, 0.8, 0.8)

        # Signal Processing

        _, tx1_settled = self.read_flag_signal(df['TX1_SETTLED'][plot_start:plot_end][::oltc_decimate])
        _, tx2_settled = self.read_flag_signal(df['TX2_SETTLED'][plot_start:plot_end][::oltc_decimate])
        _, tx3_settled = self.read_flag_signal(df['TX3_SETTLED'][plot_start:plot_end][::oltc_decimate])
        _, tx4_settled = self.read_flag_signal(df['TX4_SETTLED'][plot_start:plot_end][::oltc_decimate])

        traces = [
            ('TX Settled', tapping_flag_color, [tx1_settled, tx2_settled, tx3_settled, tx4_settled]),
        ]

        # Centerline and Label Positions
        label_positions = []
        y_positions = []

        for y in range(1, 5):
            # y_positions.append(y * major_spacing - group_offset)
            y_positions.append(y * major_spacing)
            # y_positions.append(y * major_spacing + group_offset)
            label_positions.append(y * major_spacing)

        # y_positions.append(5*major_spacing - group_offset)
        # label_positions.append(5*major_spacing - group_offset)

        ax_tx.hlines(y_positions, 0, plot_duration, linestyle='--', lw=0.6, color=center_line_color, zorder=1)

        # Plot Traces
        for (label, flag_color, flags) in traces:

            ax_tx.plot([], [], c=flag_color, label=label)

            for index, flag in enumerate(flags):
                if flag is not None:
                    for flag_segment in flag:
                        self.plot_flag(ax_tx, plot_start, flag_segment, y_positions[index], flag_color)

        self.configure_subplot(ax_tx, 'OLTC Tapping Flags', plot_duration,
                          fixed_y_lims=(group_offset, 5 * major_spacing - group_offset))

        ax_tx.set_yticks(y_positions, ['TX4', 'TX3', 'TX2', 'TX1'])
        ax_tx.grid(axis='y')
        ax_tx.tick_params(axis='y', labelsize=6.0, left=False)

        fig.align_labels()

        plt.savefig(savefile_path, bbox_inches='tight', dpi=150, format='pdf')
        # plt.savefig(savefile_path + '.pdf', bbox_inches='tight', dpi=150, format='pdf')
        # plt.savefig(savefile_path + '.png', bbox_inches='tight', dpi=150, format='png')
        plt.cla()  # Clear the current axes.
        plt.clf()  # Clear the current figure.
        plt.close()  # Closes all the figure windows.


    # def s52513_plot(self, title, load_case, test_id, scenario, nem_is_df, save_path):

        
    #     signal_layout = [
    #         [POC_SIGNALS, TERMINAL_SIGNALS, OTHER_SIGNALS],
    #         [ACTIVE_POWER_SIGNALS,BUS_VOLTAGE_SIGNALS_1,BUS_VOLTAGE_SIGNALS_2],
    #     ]

    #     num_pages = len(signal_layout)

    #     # Read timesteps from scenario eg. [0, 2.0, 10.0]
    #     time_steps: list(float) = [float(element) for element in str(scenario["Time_Steps"]).split(",")]

    #     # Specify the output PDF file path 
    #     with PdfPages(save_path) as pdf:

    #         # Loop through each page
    #         for page_number, page_signals in enumerate(signal_layout, 1):

    #             fig = self.mpl_setup(40,30)

    #             num_of_cols = len(page_signals)
    #             num_of_rows = len(page_signals[0])

    #             gs_main = gridspec.GridSpec(8, num_of_cols, figure=fig)

                
    #             cols = [
    #                 gridspec.GridSpecFromSubplotSpec(num_of_rows, 1, subplot_spec=gs_main[1:-1,0]),
    #                 gridspec.GridSpecFromSubplotSpec(num_of_rows, 1, subplot_spec=gs_main[1:-1,1]),
    #                 gridspec.GridSpecFromSubplotSpec(num_of_rows, 1, subplot_spec=gs_main[1:-1,2]),
    #             ]
                
    #             description_ax = fig.add_subplot(gs_main[0:1,0:3])
    #             test_name = str(save_path).split('\\')[-3] + ", " + str(save_path).split('\\')[-1]
    #             description_ax.set_title(f"\n{title}: {load_case}{test_id} Page {page_number} of {num_pages}", style = 'italic', fontsize='large', loc='left', y=1.1, pad=0, color=AXIS_COLOUR)
    #             description_ax.axis("off")


    #             # TABLE PLOTTING
    #             table_data = [
    #                 ['Filename/Int. Ref.', f"{scenario['File_Name']}"],
    #                 ['Test Description', f"{scenario['Category']}"],
    #             ]

    #             column_names = scenario.index.tolist()
    #             colWidths = [0.20, 0.4]

    #             description_ax.table(cellText=table_data, colWidths=colWidths,cellLoc='left', bbox=[0,0,0.4,1])

    #             def set_table_edge_color(ax, color):
    #                 table: Table
    #                 for child in ax.get_children():
    #                     if isinstance(child, Table):
    #                         table = child
    #                         break

    #                 cell: Cell
    #                 for cell in table.get_children():
    #                     if isinstance(cell, Cell):
    #                         # print(cell)
    #                         cell.set_edgecolor(AXIS_COLOUR)
    #                         cell.set_linewidth(0.75)
    #                         cell.get_text().set_color(AXIS_COLOUR)

    #             set_table_edge_color(description_ax, AXIS_COLOUR)


    #             # Windlab Logo
    #             self.plot_logo(fig)

                
    #             # Loop through each column
    #             for col_id, col_signals in enumerate(page_signals):

    #                 # Loop through each row
    #                 for row_id, signal_names in enumerate(col_signals):

    #                     # print(signal_names)

    #                     ax = fig.add_subplot(cols[col_id][row_id,-1])

    #                     for signal_index, signal_name in enumerate(signal_names):

    #                         signal_colour = SIGNAL_COLS[signal_index]

    #                         # Grab and Isolate Signals
    #                         nem_is_signal_raw = nem_is_df[signal_name][::PSSE_DECIMATE]

    #                         nem_is_t = nem_is_signal_raw.index
    #                         nem_is_v = nem_is_signal_raw.values

    #                         # Wrap PSSE phase angle 
    #                         if "DEG" in signal_name:
    #                             nem_oos_v = self.wrap_angle(nem_oos_v)
    #                             nem_is_v = self.wrap_angle(nem_is_v)

    #                         # Plot Signals
    #                         ax.plot(nem_is_t, nem_is_v, lw=SIGNAL_LINE_WIDTH, c=signal_colour, zorder=4, label=signal_name)

    #                     # Apply custom scale on V plots
    #                     if "x_axis_limits" in column_names:
    #                         xlims = [
    #                             int(str(scenario["x_axis_limits"]).split(",")[0]),
    #                             int(str(scenario["x_axis_limits"]).split(",")[1]),
    #                         ]
    #                     else:
    #                         xlims = None

    #                     self.config_subplot(ax, xlims, legend_ncol=2, signal_name=signal_name)

    #                     # Skip 1%
    #                     if "1%" in scenario["Category"]:
    #                         continue

    #                     # Settling Time Plots
    #                     analysis_values = [
    #                         ("POC_P_MW","POC_P_Settled_Time_S1","POC_P_Settling_Time_S1"),
    #                         ("POC_Q_MVAR","POC_Q_Settled_Time_S1","POC_Q_Settling_Time_S1"),
    #                         ("POC_V_PU","POC_V_Settled_Time_S1","POC_V_Settling_Time_S1"),
    #                         ("POC_P_MW","POC_P_Settled_Time_S2","POC_P_Settling_Time_S2"),
    #                         ("POC_Q_MVAR","POC_Q_Settled_Time_S2","POC_Q_Settling_Time_S2"),
    #                         ("POC_V_PU","POC_V_Settled_Time_S2","POC_V_Settling_Time_S2"),
    #                     ]

    #                     # for signal_key, settled_time_key, settling_time_key in analysis_values:
    #                     #
    #                     #     if signal_key in signal_names:
    #                     #         scenario_col_names = scenario.index.tolist()
    #                     #
    #                     #         if not settling_time_key in scenario_col_names:
    #                     #             continue
    #                     #
    #                     #         if not scenario[settling_time_key]:
    #                     #             continue
    #                     #
    #                     #         settled_time = float(scenario[settled_time_key])
    #                     #         settling_time = float(scenario[settling_time_key])
    #                     #         settled_value = float(scenario[settled_time_key.replace("Settled_Time","Settled_Value")])
    #                     #
    #                     #         ymin, ymax = ax.get_ylim()
    #                     #         yrange = ymax-ymin
    #                     #
    #                     #         ax.vlines(
    #                     #             x=settled_time,
    #                     #             ymin=settled_value + yrange*0.05 ,
    #                     #             ymax=settled_value - yrange*0.05 ,
    #                     #             colors=COL_TRIP,
    #                     #             linestyles="dashed",
    #                     #             lw = 1,
    #                     #         )
    #                     #
    #                     #         ax.text(x=settled_time+0.5, y=settled_value+0.0*yrange, s=f"Settling Time = {settling_time:.2f} s", fontsize=8)
    #                     #         ax.set_ylim(ymin=ymin-0.1*yrange,ymax=ymax+0.1*yrange)
    #                     #
    #                     #
    #                     #         # Rise Time
    #                     #         rise_time_key = settling_time_key.replace("Settling_Time","Rise_Time")
    #                     #
    #                     #         if not rise_time_key in scenario_col_names:
    #                     #             continue
    #                     #
    #                     #         if not scenario[rise_time_key]:
    #                     #             continue
    #                     #
    #                     #         rise_time = float(scenario[rise_time_key])
    #                     #
    #                     #         ax.text(x=settled_time+0.5, y=settled_value-0.08*yrange, s=f"Rise Time = {rise_time:.2f} s", fontsize=8)
    #                     #
    #                     #
    #                     # # Plot vertical lines at events for signals with analysis
    #                     # signals_with_analysis = [
    #                     #     "POC_P_MW",
    #                     #     "POC_Q_MVAR",
    #                     #     "POC_V_PU",
    #                     # ]
    #                     #
    #                     # if signal_name in signals_with_analysis:
    #                     #
    #                     #     for index, time_step in enumerate(time_steps[1:-1],start=1):
    #                     #
    #                     #         ymin, ymax = ax.get_ylim()
    #                     #         yrange = ymax-ymin
    #                     #         ax.vlines(
    #                     #             x=time_step,
    #                     #             ymin=ymin,
    #                     #             ymax=ymax,
    #                     #             colors=COL_TRIP,
    #                     #             linestyles="dashed",
    #                     #             lw = 1,
    #                     #         )
    #                     #
    #                     #         ax.text(x=time_step+0.5, y=ymax-0.08*yrange, s=f"Step  {index}", fontsize=8)
                                
                        

    #                     # Plot Title
    #                     plot_title = signal_name.replace("_"," ").replace("BESS1","BESS and WTG")
    #                     ax.set_title(plot_title, style = 'italic', fontsize='x-small', loc='left', y=1.00, pad=5, color=AXIS_COLOUR)
        
    #             pdf.savefig()   
    #             plt.close()    

    #         plt.close()




    def s5255_plot(self, title, load_case, test_id, scenario, nem_is_df, save_path):

        

        signal_layout = [
            # [POC_SIGNALS, TERMINAL_SIGNALS, OTHER_SIGNALS,],
            [ACTIVE_POWER_SIGNALS,BUS_VOLTAGE_SIGNALS_1,BUS_VOLTAGE_SIGNALS_2],
        ]

        num_pages = len(signal_layout)
        

        # Specify the output PDF file path 
        with PdfPages(save_path) as pdf:

            # Loop through each page
            for page_number, page_signals in enumerate(signal_layout, 1):

                fig = self.mpl_setup(40,30)

                num_of_cols = len(page_signals)
                num_of_rows = len(page_signals[0])

                gs_main = gridspec.GridSpec(8, num_of_cols, figure=fig)

                
                cols = [
                    gridspec.GridSpecFromSubplotSpec(num_of_rows, 1, subplot_spec=gs_main[1:-1,0]),
                    gridspec.GridSpecFromSubplotSpec(num_of_rows, 1, subplot_spec=gs_main[1:-1,1]),
                    gridspec.GridSpecFromSubplotSpec(num_of_rows, 1, subplot_spec=gs_main[1:-1,2]),
                ]
                
                description_ax = fig.add_subplot(gs_main[0:1,0:3])
                test_name = str(save_path).split('\\')[-3] + ", " + str(save_path).split('\\')[-1]
                description_ax.set_title(f"\n{title}: {load_case}{test_id} Page {page_number} of {num_pages}", style = 'italic', fontsize='large', loc='left', y=1.1, pad=0, color=AXIS_COLOUR)
                description_ax.axis("off")

                # TABLE PLOTTING

                if scenario['Fault_Bus'] == scenario['From_Bus']:
                    fault_loc = "From Bus"
                else:
                    fault_loc = "To Bus"

                fault_type_table = {
                    1: '1PHG',
                    4: '2PHG',
                    7: '3PHG',
                }

                fault_type = fault_type_table[int(scenario['Fault_Type'])]

                # TABLE PLOTTING
                table_filename = scenario['File_Name']
                table_filename = table_filename.replace(".12",".5")
                table_data = [
                    # ['Filename/Int. Ref.', f"{scenario['File_Name']}"],
                    ['Filename/Int. Ref.', table_filename],
                    ['Contingency Type', f"{scenario['Contingency_Type']}"],
                    ['Contingency Loc', f"{scenario['Contingency_Loc']}"],
                    ['Fault Type', f"{fault_type}"],
                    ['Fault End', f"{fault_loc}"],
                    ['Auto Reclose', f"{scenario['Reclose']}"],
                ]

                column_names = scenario.index.tolist()


                if "Contingency Type" in column_names:
                    
                    table_data.append(['Contingency Type', f"{scenario['Contingency_Type']}"])
                    table_data.append(['Contingency Loc', f"{scenario['Contingency_Loc']}"])

                    fault_type = fault_type_table[int(scenario['Fault_Type'])]
                    table_data.append(['Fault Type', f"{fault_type}"])

                    if scenario['Fault_Bus'] == scenario['From_Bus']:
                        fault_loc = "From Bus"
                    else:
                        fault_loc = "To Bus"

                    table_data.append(['Fault End', f"{fault_loc}"])
                    table_data.append(['Auto Reclose', f"{scenario['Reclose']}"])

                

                colWidths = [0.20, 0.4]

                description_ax.table(cellText=table_data, colWidths=colWidths,cellLoc='left', bbox=[0,0,0.3,1])

                def set_table_edge_color(ax, color):
                    table: Table
                    for child in ax.get_children():
                        if isinstance(child, Table):
                            table = child
                            break

                    cell: Cell
                    for cell in table.get_children():
                        if isinstance(cell, Cell):
                            # print(cell)
                            cell.set_edgecolor(AXIS_COLOUR)
                            cell.set_linewidth(0.75)
                            cell.get_text().set_color(AXIS_COLOUR)

                set_table_edge_color(description_ax, AXIS_COLOUR)


                # Windlab Logo
                self.plot_logo(fig)

                
                # Loop through each column
                for col_id, col_signals in enumerate(page_signals):

                    # Loop through each row
                    for row_id, signal_names in enumerate(col_signals):

                        ax = fig.add_subplot(cols[col_id][row_id,-1])

                        for signal_index, signal_name in enumerate(signal_names):
                            
                            signal_colour = SIGNAL_COLS[signal_index]

                            # Grab and Isolate Signals
                            nem_is_signal_raw = nem_is_df[signal_name][::PSSE_DECIMATE]

                            nem_is_t = nem_is_signal_raw.index
                            nem_is_v = nem_is_signal_raw.values

                            # Wrap PSSE phase angle 
                            if "DEG" in signal_name:
                                nem_oos_v = self.wrap_angle(nem_oos_v)
                                nem_is_v = self.wrap_angle(nem_is_v)

                            # right axis for P
                            if signal_name in ["WT1_P_MW","WT1_Q_MVAR"]:
                                # Create a twinx() to add a secondary y-axis on the right
                                ax2 = ax.twinx()
                                line2 = ax2.plot(nem_is_t, nem_is_v, lw=SIGNAL_LINE_WIDTH, c=signal_colour, zorder=4, label=signal_name)
                                # ax2.set_ylabel('Right Y-axis', color='r')
                                # ax2.tick_params(axis='y', labelcolor='r')
                                ax2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
                                ax2.tick_params(axis='y', labelsize=8)
                                ax2.set_ylim(min(-2, round(min(nem_is_v)*1.35)), max(8,round(max(nem_is_v)*1.35)))

                                ax2.set_ylabel(f"{signal_name}",color=line2[0].get_color(), fontsize='small')
                                # ax2.set_ylabel(f"{signal_name}", fontsize='x-small')
                                ax2.yaxis.set_label_coords(0.97,0.5)
                                ax2.legend(frameon=False, fontsize='x-small', bbox_to_anchor=(0.8, 1.1),
                                          ncol=1,
                                          borderpad=0, columnspacing=1, handletextpad=0.3, handlelength=1.2,
                                          labelcolor=AXIS_COLOUR)

                            elif signal_name in ["BESS1_P_MW","BESS1_Q_MVAR"]:
                                line1 = ax.plot(nem_is_t, nem_is_v, lw=SIGNAL_LINE_WIDTH, c=signal_colour, zorder=4,
                                        label=signal_name)
                                ax.set_ylabel(f"{signal_name}", color=line1[0].get_color(), fontsize='small')
                                ax.yaxis.set_label_coords(0.03, 0.5)

                            elif signal_name =="FPOC_HZ":
                                ax2 = ax.twinx()
                                line2 = ax2.plot(nem_is_t, nem_is_v, lw=SIGNAL_LINE_WIDTH, c=signal_colour, zorder=4,
                                                 label=signal_name)
                                ax2.legend(frameon=False, fontsize='x-small', bbox_to_anchor=(1.0, 1.1),
                                           ncol=1,
                                           borderpad=0, columnspacing=1, handletextpad=0.3, handlelength=1.2,
                                           labelcolor=AXIS_COLOUR)

                                # Remove the left y-axis labels and ticks
                                ax.set_yticks([])

                            else:
                            # Plot Signals
                                ax.plot(nem_is_t, nem_is_v, lw=SIGNAL_LINE_WIDTH, c=signal_colour, zorder=4, label=signal_name)

                        # Apply custom scale on V plots
                        if "x_axis_limits" in column_names:
                            xlims = [
                                int(str(scenario["x_axis_limits"]).split(",")[0]),
                                int(str(scenario["x_axis_limits"]).split(",")[1]),
                            ]
                        else:
                            xlims = None
                        

                        if signal_name in ["WT1_P_MW","WT1_Q_MVAR","FPOC_HZ"]:
                            self.config_subplot(ax2, xlims, legend_ncol=2, signal_name=signal_name)
                        else:
                            self.config_subplot(ax, xlims, legend_ncol=2, signal_name=signal_name)


                        # Active Power Recovery
                        if signal_name == "POC_P_MW":
                            # if scenario["Ppoc_95pc_Recovery_Time"] == "":
                            #     break
                            #
                            # recovery_duration = float(scenario["Ppoc_95pc_Recovery_Time"])
                            # fault_time = float(scenario["Fault_Time"])
                            # fault_duration = float(scenario["Fault_Duration"])
                            # recovery_time = fault_time + fault_duration + recovery_duration

                            ymin, ymax = ax.get_ylim()
                            if abs(ymax-ymin) <200:
                                ymin = ymin-80
                                ymax = ymax+80
                            # ax.vlines(
                            #     x=fault_time,
                            #     ymin=ymin,
                            #     ymax=ymax,
                            #     colors=COL_TRIP,
                            #     linestyles="dashed",
                            #     lw = 1,
                            # )
                            # ax.vlines(
                            #     x=recovery_time,
                            #     ymin=ymin,
                            #     ymax=ymax,
                            #     colors=COL_TRIP,
                            #     linestyles="dashed",
                            #     lw = 1,
                            # )
                            # ax.text(x=recovery_time+0.5, y=ymin+20, s=f"95% Recovery Time = {recovery_duration:.2f} s", fontsize=8)
                            ax.set_ylim(ymin=ymin,ymax=ymax)

                        # scale plot 77 POC V
                        if (signal_name == "POC_V_PU"):
                            ymin, ymax = ax.get_ylim()
                            if abs(ymax-ymin) < 0.01:
                                ax.set_ylim(ymin=ymin-0.3, ymax=ymax+0.3)

                        # scale plot 77 POC Q
                        if (signal_name == "POC_Q_MVAR"):
                            ymin, ymax = ax.get_ylim()
                            if abs(ymax - ymin) < 20:
                                ax.set_ylim(ymin=ymin - 20, ymax=ymax + 20)

                        # scale plot 77 BESS1 and WTG1 Q
                        if (signal_name == "BESS1_Q_MVAR"):
                            ymin, ymax = ax.get_ylim()
                            if abs(ymax - ymin) < 20:
                                ax.set_ylim(ymin=ymin - 20, ymax=ymax + 20)

                        # scale plot 77 BESS1 and WTG1 V
                        if (signal_name == "BESS1_V_PU"):
                            ymin, ymax = ax.get_ylim()
                            if abs(ymax - ymin) < 0.01:
                                ax.set_ylim(ymin=ymin - 0.3, ymax=ymax + 0.3)

                        # scale plot 77 BESS1 and WTG1 IQ
                        if (signal_name == "BESS1_IQ_PU"):
                            ymin, ymax = ax.get_ylim()
                            if abs(ymax - ymin) < 0.01:
                                ax.set_ylim(ymin=ymin - 0.3, ymax=ymax + 0.3)

                        # Plot Title
                        plot_title = signal_name.replace("_"," ").replace("BESS1","BESS and WTG")
                        ax.set_title(plot_title, style = 'italic', fontsize='x-small', loc='left', y=1.00, pad=5, color=AXIS_COLOUR)
        
                pdf.savefig()   
                plt.close()    

            plt.close()


    def find_pickle(self, dir, key_str):
        file_list = []

        if not os.path.exists(dir):
            if DEBUG: print(f"Cannot locate: \n {dir}")
            return None

        for file_name in os.listdir(dir):
            # print(file_name)
            if file_name.endswith(".pkl") and key_str in file_name:
                file_list.append(file_name)

        if len(file_list) > 1 or  len(file_list) > 1:
            if DEBUG: print("More than 1 matching file found.")
            return None

        if len(file_list) == 0:
            if DEBUG: print(f"{key_str} not found in {dir}") 
            return None

        pickle_file = file_list[-1]
        # print(pickle_file)
    
        return pickle_file
    

    def _get_filenames(self, dir):
        filenames = []
        for file_name in os.listdir(dir):
            if os.path.isfile(os.path.join(dir, file_name)):
                filename_without_ext = os.path.splitext(file_name)[0]
                filenames.append(filename_without_ext)
        return filenames
    

    def mpl_setup(self, w, h):

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
        plt.subplots_adjust(left=0.04, right=0.975, bottom=0.05, top=0.92, wspace=0.1, hspace=0.3)

        return fig


    def wrap_angle(self, signal):

        new_signal = signal

        for i, value in enumerate(signal):

            angle = value % 360
            angle = (angle + 360) % 360

            if (angle > 180):
                angle -= 360

            new_signal[i] = angle

        return new_signal
    
    def plot_logo(self, fig):
            gs_logo = gridspec.GridSpec(12, 8, figure=fig)
            logo_ax = fig.add_subplot(gs_logo[0,-1])
            logo = plt.imread("./windlab.jpg")
            logo_ax.imshow(logo)
            logo_ax.axis("off")
    
    def config_subplot(self, ax,xlims,legend_y=1.10, legend_ncol=2, signal_name=""):

        initial_xlims = ax.get_xlim()
        initial_ylims = ax.get_ylim()
        

        # Grid Properties:
        ax.set_axisbelow(True)
        ax.grid(which='major', linestyle='--', linewidth=0.5, color='black', alpha=0.4, zorder=-100)
        ax.grid(which='minor', linestyle='--', linewidth=0.5, color='black', alpha=0.1, zorder=-100)

        #  Creates a XMajor Loactor that attempts to place
        xmajor_locator = matplotlib.ticker.MaxNLocator(nbins=5, steps=[1, 2, 4, 5, 10])
        
        if xlims:
            ax.set_xlim(xlims[0], xlims[1])
            
            major_xticks_values = xmajor_locator.tick_values(xlims[0], xlims[1])
        else:       
            major_xticks_values = xmajor_locator.tick_values(9, initial_xlims[1])
       
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


        current_y_min, current_y_max = initial_ylims
        current_y_range = current_y_max - current_y_min
        min_y_range = 0.05

        if "_V_PU" in signal_name:
            if signal_name not in [
                "GMSS_275KV_V_PU", 
                "ROSS_275KV_V_PU",
                "KANBAN_POC_275KV_V_PU",
                "CHALUMBIN_275KV_V_PU",
                "STRATHMORE_275KV_V_PU",
                "YABULLU_STH_275KV_V_PU",
                "TULLY_275KV_V_PU",
                "KIDPUMHYD_G1_16.5KV_V_PU",
            ]:
                
                min_y_range = 0.00001 # disable it
                ax.autoscale(enable=True, axis="y", tight=True)
                ax.margins(0.02, 0.02)
                
            else:
                min_y_range = 0.04
                
            

        if "_P_MW" in signal_name:
            min_y_range = 5

        if "_HZ" in signal_name:
            min_y_range = 1

        if "_DEG" in signal_name:
            min_y_range = 10


        if current_y_range < min_y_range:

            y_min = ((current_y_min + current_y_max) - min_y_range)/2
            y_max = ((current_y_min + current_y_max) + min_y_range)/2
            ax.set_ylim(y_min,y_max)

        # Do not limit x axis for interconnector signals
        if signal_name in [
            "QNI1_P_MW", 
            "QNI2_P_MW",
            "TERRANORA1_P_MW",
            "TERRANORA2_P_MW",
        ]:  
            ax.set_xlim(9, initial_xlims[1])

        # Y-axis    
        ax.tick_params(axis='y',  colors=AXIS_COLOUR, labelsize=8, direction='in', which='both')

        # Legend.
        # if show_legend:  # bbox to anchor 1.05
        ax.legend(frameon=False, fontsize='x-small', bbox_to_anchor=(1.0, legend_y), ncol=legend_ncol,
                borderpad=0, columnspacing=1, handletextpad=0.3, handlelength=1.2, labelcolor=AXIS_COLOUR)
            
        ax.spines["bottom"].set_color(AXIS_COLOUR)      
        ax.spines["top"].set_color(AXIS_COLOUR)
        ax.spines["left"].set_color(AXIS_COLOUR)
        ax.spines["right"].set_color(AXIS_COLOUR)

    def configure_subplot(self, ax, subplot_title, plot_duration, legendon=True, xlabelson=False, min_y_range=None,
                          fixed_y_lims=None):
        # Grid Properties:
        ax.set_axisbelow(True)
        ax.grid(which='major', linestyle='--', linewidth=0.5, color='darkgray', alpha=0.8, zorder=-100)
        ax.grid(which='minor', linestyle='--', linewidth=0.5, color='darkgray', alpha=0.4, zorder=-100)

        # X-Axis Propertie:
        ax.set_xlim(0, plot_duration)
        #  Creates a XMajor Loactor that attempts to place
        xmajor_locator = matplotlib.ticker.MaxNLocator(nbins=5, steps=[1, 2, 4, 5, 10])
        major_xticks_values = xmajor_locator.tick_values(0, plot_duration)
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

        # Legend.
        if legendon:
            ax.legend(frameon=False, fontsize='x-small', bbox_to_anchor=(1.0, 0.0), ncol=6,
                      borderpad=0, columnspacing=1, handletextpad=0.3, handlelength=1.2, labelcolor=AXIS_COLOUR)

        ax.spines["bottom"].set_color(AXIS_COLOUR)
        ax.spines["top"].set_color(AXIS_COLOUR)
        ax.spines["left"].set_color(AXIS_COLOUR)
        ax.spines["right"].set_color(AXIS_COLOUR)

        ax.set_xlim(0, plot_duration)

    def signal_plot(self, ax, title, traces, start, duration, fixed_y_lims=None, min_y_range=None, time_axis_on=False):
        if time_axis_on:
            legendon = False
            xlabelson = True
        else:
            legendon = True
            xlabelson = False

        for signal_label, lw, colour, order, signal in traces:
            ax.plot(signal.index - start, signal.values,
                    lw=lw, c=colour, zorder=order, label=signal_label)

        if fixed_y_lims is not None:
            self.configure_subplot(ax, title, duration, fixed_y_lims=fixed_y_lims, legendon=legendon, xlabelson=xlabelson)
        elif min_y_range is not None:
            self.configure_subplot(ax, title, duration, min_y_range=min_y_range, legendon=legendon, xlabelson=xlabelson)
        else:
            self.configure_subplot(ax, title, duration, legendon=legendon, xlabelson=xlabelson)

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

    def read_psse_frt_flags(self, flag_signal):
        """
        Description:

        Function to read and interpret goldwind ppc psse flag signals, returning
        a list of intervals at which the flag is high or low

        Inputs: Flag signal that is either:

            0 Normal, 1 Prefault, 2 LVRT, 3 Recoverage, -2 HVRT"

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

            # difference between current and previous flag value
            difference = value - previous_value

            # start LVRT
            if previous_value != value and value == 2:
                low_flag_start = time

            # stop LVRT
            if previous_value == 2 and value != 2:
                low_flag_stop = time
                low_flag.append(pd.Series(data=[1.0, 1.0], index=[low_flag_start, low_flag_stop]))

            # start HVRT
            if previous_value != value and value == -2:
                high_flag_start = time

            # stop HVRT
            if previous_value == -2 and value != -2:
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
            # marker = '',
            # markersize=0.8,
            lw=5,
            solid_capstyle='butt',
            c=flag_color,
            zorder=2,
        )

    def low_pass_filter(self, signal: pd.Series, cut_off: float):
        fs = 1 / (signal.index[1] - signal.index[0])
        b, a = butter(1, cut_off, fs=fs, btype='low', analog=False)
        filtered_signal = lfilter(b, a, signal.values - signal.values[0]) + signal.values[0]

        return pd.Series(
            data=filtered_signal,
            index=signal.index,
        )


if __name__ == "__main__":

    bmu = NEMPostprocessUtility()
