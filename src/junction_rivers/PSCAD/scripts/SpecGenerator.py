import pandas as pd
import json
from typing import Dict
import math
from icecream import ic
import random

class CalcSheetError(Exception):
    def __init__(self):
        print(f"Calc Sheet Error")
    
    def __init__(self, input):
        print(f"Calc Sheet Error: {input}")
class SpecGenerator():
    
    HEADERS = [
        'Grouping',
        'Key_Tests',
        'File_Name',
        'Post_Init_Duration_s',
        'Category',
        'DIR',
        ]
    
    WTG_PZERO = 0
    BESS_PMAX = 1
    BESS_PZERO = 2
    BESS_PMIN = 3
    
    def __init__(self, calc_sheet_path, spec_path):
        self.spec_generator(calc_sheet_path, spec_path)
        ic(calc_sheet_path)
        ic(spec_path)
        
                
    def spec_generator (self, calc_sheet_path, spec_path):
        
        checklist = pd.read_excel(calc_sheet_path, sheet_name="Checklist", index_col="Checklist", header=0)
        ic(checklist)
        self.system_inf = pd.read_excel(calc_sheet_path, sheet_name="System_inf", index_col="Var name", header=0)
        ic(self.system_inf)
        
        # Record project information
        self.v_nom = self.system_inf.loc["Nominal Voltage"]["Var Val"]
        self.f_nom = self.system_inf.loc["Nominal Frequency"]["Var Val"]
        self.p_nom = self.system_inf.loc["Nominal Power"]["Var Val"]
        self.q_nom = self.p_nom * 0.395
        self.s_nom = math.sqrt(self.p_nom**2 + self.q_nom**2)
        self.z_base = self.v_nom**2/self.p_nom
        self.q_set = 0
        self.v_set = 1.034
        self.p_max_bess = 70

        # Initialise data frame for the output spec data
        spec_df = pd.DataFrame()
        
        # Iterrate through the categories as recorded in the checklist
        for index_cat, category in checklist.iterrows():
            ic(index_cat)
            ic(category)
        
            if category.loc["Action"] == "Yes":
                tests = pd.read_excel(calc_sheet_path, sheet_name=index_cat, header=0)
                ic(tests)
                
                # Iterate through the tests in the sheet
                for _, row in tests.iterrows():
                    ic(row)
                    
                    if row["Action"] =="Yes":
                        category_formated = index_cat.lower().replace(" ", "_")
                        new_row = dict()
                        new_row.update({'File_Name': f'{category_formated}_test_{row["Test"]}'.replace("-","_")})
                        new_row.update({'DIR': category_formated})
                        new_row.update({'En_SMIB_init_v': 1})
                        new_row.update({'Infinite_Grid_v': 1})
                        new_row.update(self.add_new_row(row))
                        ic(new_row)
                        spec_df = spec_df.append(new_row, ignore_index=True)
                        ic(spec_df)
                        
        spec_df.to_csv(spec_path)
        
        
################# ADD PARAMS FOR A TEST ################

    def add_new_row(self, row: pd.DataFrame):
        ic('add_new_row')
        new_row=dict()
        ic(new_row)
        new_row.update(self.add_file_info(row))
        ic(new_row)
        new_row.update(self.add_v_specs(row))
        ic(new_row)
        new_row.update(self.add_q_specs(row))
        ic(new_row)
        new_row.update(self.add_p_specs(row, self.BESS_PZERO))
        ic(new_row)
        new_row.update(self.add_freq_specs(row))
        ic(new_row)
        new_row.update(self.add_fault_specs(row))
        
        return new_row
    
    def add_file_info(self, row: pd.DataFrame):
        ic('add_file_info')
        row_sect = dict()
        row_sect =  {
                    'Grouping': '',
                    'Key_Tests': '',
                    'Category': '',
                    'Post_Init_Duration_s': row["End Run (s)"],
                    'Grid_SCR_v': self.read_scr(row),
            }
        ic(row_sect)
        return row_sect
    
    def add_v_specs(self, row: pd.DataFrame):
        ic('add_v_specs')
        row_sect = dict()
        if 'Vgrid' in row:
            ic("Vgrid in row")
            row_sect = {'Init_Vpoc_pu_v': self.calc_vpoc_from_vslack(row),
                        'Vslack_pu_v': row["Vgrid"]}
        elif 'Vpoc' in row:
            ic("Vpoc in row")
            row_sect = {'Init_Vpoc_pu_v': row["Vpoc"],
                        'Vslack_pu_v': self.calc_vslack_from_vpoc(row)}
        else:
            raise CalcSheetError
        ic(row_sect)
        return row_sect
    
    def add_q_specs(self, row: pd.DataFrame):
        ic('add_q_specs')
        row_sect = dict()
        row_sect = {'Init_Qpoc_pu_v': self.q_set}
        ic(row_sect)
        return row_sect
    
    def add_p_specs(self, row: pd.DataFrame, wf_state: int):
        ic('add_p_specs')
        row_sect = dict()
        if 'Active Power (pu)' in row:
            if wf_state == self.WTG_PZERO:
                row_sect = {'Pref_Wind_MW_v': 0,
                        'Pref_BESS_MW_v': row["Active Power (pu)"]*self.p_nom}
            elif wf_state == self.BESS_PMAX:
                row_sect = {'Pref_Wind_MW_v': row["Active Power (pu)"]*self.p_nom,
                        'Pref_BESS_MW_v': self.p_max_bess}
            elif wf_state == self.BESS_PZERO:
                row_sect = {'Pref_Wind_MW_v': row["Active Power (pu)"]*self.p_nom,
                        'Pref_BESS_MW_v': 0}
            elif wf_state == self.BESS_PMIN:
                row_sect = {'Pref_Wind_MW_v': row["Active Power (pu)"]*self.p_nom,
                        'Pref_BESS_MW_v': -self.p_max_bess}
            else:
                row_sect = {'Pref_Wind_MW_v': row["Active Power (pu)"]*self.p_nom/2,
                        'Pref_BESS_MW_v': row["Active Power (pu)"]*self.p_nom/2}
        else:
            row_sect = {'Pref_Wind_MW_v': 0,
                        'Pref_BESS_MW_v': 0}
        ic(row_sect)
        return row_sect
    
    def add_freq_specs(self, row: pd.DataFrame):
        ic('add_feq_specs')
        row_sect = dict()
        if 'Freq target' in row:
            if 'Apply step (s)' in row:
                if 'Freq ramp Hz/s':
                    row_sect = {'Grid_Hz_t': f'[0, {row["Apply step (s)"]}, {self.calc_ramp(row)}, {row["End Run (s)"]}]',
                                'Grid_Hz_v': f'[{self.f_nom}, {row["Freq target"]}, {row["Freq target"]}]',
                                'Grid_Hz_r': f'[true, false]'}
                else:
                    row_sect = {'Grid_Hz_t': f'[0, {row["Apply step (s)"]}, {row["End Run (s)"]}]',
                                'Grid_Hz_v': f'[{self.f_nom}, {row["Freq target"]}]'}
        else:
            row_sect = {'Grid_Hz_v': f'[{self.f_nom}]'}
        ic(row_sect)
        return row_sect
    
    def add_fault_specs(self, row: pd.DataFrame):
        ic("add_fault_specs")
        row_sect = dict()
        if 'Apply Fault (s)' in row:
            if 'Fault Impedance' in row:
                if 'Fault duration (s)' in row:
                    if row["Fault duration (s)"].isnumeric():
                        if 'Fault type' in row:
                            row_sect = {'Init_Fault_MVA': self.calc_fault_level(row),
                                        'Init_Fault_X_on_R': self.read_x2r(row),
                                        'Fault_Type_v': self.get_fault_code(row["Fault type"]),
                                        'Fault_Depth': self.calc_fault_voltage(row),
                                        'Fault_Time': row["Apply Fault (s)"],
                                        'Fault_Duration': row["Fault duration (s)"],
                                        'Post_Fault_Duration': row["End Run (s)"] - row["Apply Fault (s)"] - row["Fault duration (s)"],
                                        'Enable_Fault_Studies_v': 1,
                                        'Fault_Strategy_v': 0,
                                        'Fault_X2R_v': self.read_x2r(row),
                                        'Fault_Timing_Signal_t': f'[0, {row["Apply Fault (s)"]}, {row["Apply Fault (s)"] + row["Fault duration (s)"]}',
                                        'Fault_Timing_Signal_v': f'[0, 1, 0]',
                                        'Ures_v': self.calc_fault_voltage(row)}
                    elif row["Fault duration (s)"] in ["S1", "S2", "S3", "S4", "S5"]:
                        row_sect = {'Init_Fault_MVA': self.calc_fault_level(row),
                                    'Init_Fault_X_on_R': self.read_x2r(row),
                                    'Fault_X2R_v': self.read_x2r(row),
                                    'Enable_Fault_Studies_v': 1,
                                    'Fault_Strategy_v': 0,}
                        row_sect.update(self.mfrt_seq(row))
        ic(row_sect)                    
        return row_sect
        
##################### MORE TOOLS #######################

    def mfrt_seq(self, row: pd.DataFrame):
        ic('mfrt_seq')
        random.seed(1)
        fault_info = []
        fault_type_options = ["3PHG", "2PHG", "1PHG", "L-L"]
        if row["Fault duration (s)"] == "S1":
            fault_info = [("3PHG", 5, 0.1, 0.25),
                        ("3PHG", 5.25, 0.1, 0.25),
                        ("3PHG", 5.5, 0.1, 0.25),
                        ("2PHG", 8, 0.1, 0.25),
                        ("2PHG", 11, 0.1, 0.25),
                        ("2PHG", 13, 0.1, 0.25)]
            time = [13]
            fault_duration = 0.1
        else:
            # From S5.2.5.5 consider up to 15 faults in up to 5 mins.
            # Arbitrarily I decided to constrain it to between 5 and 15 faults in 2 to 5 mins.
            min_no_faults = 5
            max_no_faults = 15
            max_time = 300
            # Choose a ranom number of faults to apply between 5 and 15
            no_faults = random.randint(min_no_faults, max_no_faults)
            ic(no_faults)
            # We can specify a minimum time between faults
            min_time_between_faults = 0.1
            # By using the integer random sampling method we can control the accuracy.
            accuracy = 0.001
            time_range = range(0, int(min(max_time, row["End Run (s)"])/accuracy), int(min_time_between_faults/accuracy))
            ic(time_range)
            time = random.sample(time_range, k=no_faults)
            time.sort()
            # Convert the list back to seconds.
            time = [item*accuracy for item in time]
            ic(time)
            # Loop through each time step.
            for i in range(len(time)):
                fault_type = random.choice(fault_type_options)
                next_time = int(row["End Run (s)"]) if i == len(time)-1 else time[i+1]
                fault_duration = random.randrange(1, int((next_time-time[i])/accuracy), 1)*accuracy
                fault_multiplier = random.randint(0, 5)
                fault_info.append((fault_type, time[i], fault_duration, fault_multiplier))
        row_sect = self.generate_mfrt(fault_info)
        ic(row_sect)
        row_sect['Post_Fault_Durations'].extend([row["End Run (s)"] - time[-1] - fault_duration])
        ic(row_sect)
        return row_sect
    
    def generate_mfrt(self, fault_info: list):
        ic('generate_mfrt')
        row_sect = dict()
        fault_types = []
        fault_durations = []
        post_fault_durations = []
        fault_zs_multiplier = []
        fault_timing_signal_v = [0]
        fault_timing_signal_t = [0]
        fault_times = []
        
        for index, (fault_type, fault_time, fault_duration, fault_impedence_multiplier) in enumerate(fault_info):

            fault_types.append(self.get_fault_code(fault_type))
            fault_durations.append(fault_duration)
            ic(fault_durations)
            if index > 0 and index < len(fault_info) - 1:
                post_fault_durations.append(fault_time - fault_times[index-1] - fault_durations[index-1])
            ic(post_fault_durations)
            fault_zs_multiplier.append(fault_impedence_multiplier)
            fault_timing_signal_t.extend([fault_time, fault_time+fault_duration])
            fault_timing_signal_v.extend([1, 0])
            fault_times.append(fault_time)
            ic(fault_times)
        
        row_sect = {'Fault_Types': fault_types,
                    'Fault_Durations': fault_durations,
                    'Post_Fault_Durations': post_fault_durations,
                    'Fault_Zs_Multiplier': fault_zs_multiplier,
                    'Fault_Timing_Signal_t': fault_timing_signal_t,
                    'Fault_Timing_Signal_v': fault_timing_signal_v,
                    'Zf2Zs_t': fault_times,
                    'Zf2Zs_v': fault_zs_multiplier,
                    'Fault_ Type_t': fault_times,
                    'Fault_Type_v': fault_types}    
        return row_sect
    
    def calc_ramp(self, row: pd.DataFrame):
        ic('calc_ramp')
        d_freq = abs(row["Freq target"] - self.f_nom)
        d_t = d_freq/row["Freq ramp Hz/s"]
        t_end_ramp = min(row["Apply step (s)"] + d_t, row["End Run (s)"])
        return t_end_ramp
        
        
##################### FAULT CALCULATIONS #######################
    def calc_fault_impedance(self, row: pd.DataFrame):
        ic('calc_fault_impedance')
        min_fault_impedance = 0.000000001
        if row["Fault Impedance"].isnumeric():
            fault_impedance = row["Fault Impedance"]
        elif row["Fault Impedance"][0:3] =="Zf=":
            if row["Fault Impedance"] == "Zf=0":
                fault_impedance = min_fault_impedance
            elif row["Fault Impedance"] == "Zf=Zs":
                fault_impedance = self.calc_fault_impedance(row, 1)
            else:
                raise CalcSheetError
        return fault_impedance
    
    # def calc_fault_multiplier(self, string: str):
    #     multiplier = 0
    #     if string == "Zf=0":
    #         multiplier = 0
    #     elif string == "Zf=Zs":
    #         multiplier = 1
    #     return multiplier
    
    def calc_fault_voltage(self, row: pd.DataFrame):
        ic('calc_fault_voltage')
        Zf = self.calc_fault_impedance(row)
        (_, _, Zs) = self.calc_grid_impedence(row, 1)
        Udip = Zf/(Zs + Zf)
        return Udip

# ##
# 1: phase A to ground
# 2: phase B to ground
# 3: phase C to ground
# 4: phase A, B to ground
# 5: phase A, C to ground
# 6: phase B, C to ground
# 7: phase A, B, C to ground
# 8: phase A to B
# 9: phase A to C
# 10: phase B to C
# ##
    def get_fault_code(self, code: str):
        ic('get_fault_code')
        fault_code = 0
        if code == "3PHG":
            fault_code = 7
        elif code == "2PHG":
            fault_code = 4
        elif code == "1PHG":
            fault_code = 1
        elif code == "L-L":
            fault_code = 8
        return fault_code
        
##################### CALCULATIONS #######################
    def read_scr(self, row: pd.DataFrame):
        ic("read_scr")
        scr = 0
        if row["SCR"].isnumeric():
            scr = row["SCR"]
        elif row["SCR"] == "POC":
            scr = json.loads(self.system_inf.loc["POC SCR"]["Var Val"])[1]
        else:
            raise CalcSheetError("Cannot identify SCR entry")
        ic(scr)
        return scr
        
        
    def read_x2r(self, row: pd.DataFrame):
        ic("read_scr")
        x2r = 0
        if row["X/R"].isnumeric():
            x2r = row["X/R"]
        elif row["X/R"] == "POC":
            x2r = self.system_inf.loc["POC XR ratio"]["Var Val"]
        else:
            CalcSheetError("Cannot identify X/R entry")
        ic(x2r)
        return x2r
    
    def calc_fault_level(self, row: pd.DataFrame):
        ic("calc_fault_level")
        scr = self.read_scr(row)
        fault_level = scr * self.p_nom
        ic(fault_level)
        return fault_level
    
    def calc_grid_impedence(self, row: pd.DataFrame, pu: bool):
        ic("calc_grid_impedence")
        x2r = self.read_x2r(row)
        fl = self.calc_fault_level(row)
        z = self.v_nom**2 / fl
        r = z / math.sqrt(1+x2r**2)
        x = r * x2r
        ic("ohms")
        ic(r)
        ic(x)
        l = x / (2 * self.f_nom * math.pi)
        ic(l)
        if pu == True:
            r = r/self.z_base
            x = x/self.z_base
            z = z/self.z_base
            ic("pu")
            ic(r)
            ic(x)
        l = x / (2 * self.f_nom * math.pi)
        ic(l)
        return (r, x, z)
        
    def calc_vslack_from_vpoc(self, row: pd.DataFrame):
        ic("calc_vslack_from_vpoc")
        vpoc = row["Vpoc"]
        ppoc = row["Active Power (pu)"]
        qpoc = self.q_set
        spoc = math.sqrt(ppoc**2 + qpoc**2)
        (rgrid_pu, xgrid_pu, zgrid_pu) = self.calc_grid_impedence(row, 1)
        vslack = math.sqrt(vpoc**2 + spoc**2/vpoc**2*zgrid_pu**2 - 2*(ppoc*rgrid_pu + qpoc*xgrid_pu))
        return vslack
        
    def calc_vpoc_from_vslack(self, row: pd.DataFrame):
        ic("calc_vpoc_from_vslack")
        vslack = row["Vgrid"]
        ppoc = row["Active Power (pu)"]
        qpoc = self.q_set
        spoc = math.sqrt(ppoc**2 + qpoc**2)
        (rgrid, xgrid, zgrid) = self.calc_grid_impedence(row, 1)
        rgrid_pu = rgrid/ self.z_base
        xgrid_pu = xgrid/ self.z_base
        zgrid_pu = zgrid/ self.z_base
        vgrid = math.sqrt((vslack**2 + 2*(ppoc*rgrid_pu + qpoc*xgrid_pu) + math.sqrt((vslack**2 + 2*(ppoc*rgrid_pu + qpoc*xgrid_pu))**2 - 4*spoc**2*zgrid_pu**2))/2)
        ic(vgrid)
        return vgrid
    
    