import pandas as pd
import json
from typing import Dict
import math
from icecream import ic
import random
from copy import deepcopy

class CalcSheetError(Exception):
    def __init__(self):
        print(f"Calc Sheet Error")
    
    def __init__(self, input):
        print(f"Calc Sheet Error: Cannot interpret {input}")
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
        self.calc_sheet_path = calc_sheet_path
        self.spec_generator(calc_sheet_path, spec_path)
        ic(calc_sheet_path)
        ic(spec_path)
        
                
    def spec_generator (self, calc_sheet_path, spec_path):
        
        checklist = pd.read_excel(calc_sheet_path, sheet_name="Checklist", index_col="Checklist", header=0)
        ic(checklist)
        self.system_inf = pd.read_excel(calc_sheet_path, sheet_name="System_inf", index_col="Var name", header=0)
        ic(self.system_inf)
        
        # Record project information
        self.v_nom = self.system_inf.loc["POC base Voltage"]["Var Val"]
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
                    
                    if row["Action"] =="Yes":
                        ic(row)
                        category_formated = index_cat.lower().replace(" ", "_")
                        new_row = self.add_new_row(row, category=category_formated)
                        ic(new_row)
                        for row in new_row:
                            spec_df = spec_df.append(row, ignore_index=True)
                        ic(spec_df)
                        
        spec_df.to_csv(spec_path)
        
        
################# ADD PARAMS FOR A TEST ################
    def update_rows(self, specs: list, in_row: dict):
        ic("update_rows")
        ic(specs)
        # create a new list of tests
        out_rows = []
        # for multiple specs, we want to add new tests for each one
        for index, spec in enumerate(specs):
            # create a copy of the existing test to update without changing the original one
            ic(in_row)
            temp = deepcopy(in_row)
            ic(temp)
            ic(type(temp))
            temp.update({'File_Name': f'{in_row["File_Name"]}-{index}'})
            temp.update(spec)
            out_rows.append(temp) 
        return out_rows          
    
    def add_new_row(self, row: pd.DataFrame, category: str):
        ic('add_new_row')
        # Initialise a dict with category information and a list of dicts with just the empty dict.
        first_test = dict()
        first_test.update({'File_Name': f'{category}_test_{row["Test"]}'.replace("-","_")})
        first_test.update({'DIR': category})
        first_test.update({'En_SMIB_init_v': 1})
        first_test.update({'Infinite_Grid_v': 1})
        new_tests = [first_test]
        
        file_infos = self.add_file_info(row)
        new_tests = self.update_rows(file_infos, first_test)
        ic(new_tests)
        
        out_rows = []
        for new_test in new_tests:
            q_specs = self.add_q_specs()
            out_rows.extend(self.update_rows(q_specs, new_test))
        new_tests = out_rows[:]
        ic(new_tests)
        
        out_rows = []
        for new_test in new_tests:
            p_specs = self.add_p_specs(row, self.BESS_PZERO)
            out_rows.extend(self.update_rows(p_specs, new_test))
        new_tests = out_rows[:]
        ic(new_tests)
        
        out_rows = []
        for new_test in new_tests:
            v_specs = self.add_v_specs(row, new_test)
            out_rows.extend(self.update_rows(v_specs, new_test))
        new_tests = out_rows[:]
        ic(new_tests)
        
        out_rows = []
        for new_test in new_tests:
            freq_specs = self.add_freq_specs(row)
            out_rows.extend(self.update_rows(freq_specs, new_test))
        new_tests = out_rows[:]
        ic(new_tests)
        
        out_rows = []
        for new_test in new_tests:
            fault_specs = self.add_fault_specs(row, new_test)
            ic(fault_specs)
            out_rows.extend(self.update_rows(fault_specs, new_test))
        new_tests = out_rows[:]
        ic(new_tests)
        
        return new_tests
    
    def add_file_info(self, row: pd.DataFrame):
        ic('add_file_info')
        row_sect_list = []
        row_sect = dict()
        row_sect =  {
                    'Grouping': '',
                    'Key_Tests': '',
                    'Category': '',
                    'Post_Init_Duration_s': row["End Run (s)"],
                }
        scr_and_x2rs = self.read_scr_and_x2r(row)
        for scr_and_x2r in scr_and_x2rs:
            (scr, x2r) = scr_and_x2r
            temp_sect = row_sect.copy()
            temp_sect.update({
                'Grid_SCR_v': scr,
                'Grid_X2R_v': x2r,
                'Grid_MVA_v': self.calc_fault_level(scr)
            })
            row_sect_list.append(temp_sect)
        ic(row_sect_list)
        return row_sect_list
    
    def add_q_specs(self):
        ic('add_q_specs')
        row_sect = dict()
        row_sect = {'Init_Qpoc_pu_v': self.q_set}
        ic(row_sect)
        return [row_sect]
    
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
        return [row_sect]
    
    def add_v_specs(self, row: pd.DataFrame, new_tests: list):
        ic('add_v_specs')
        row_sect = dict()
        # Use data which has already been added to this test to calculate the values to be added
        ppoc = new_tests["Pref_Wind_MW_v"] + new_tests["Pref_BESS_MW_v"]
        qpoc = new_tests["Init_Qpoc_pu_v"]
        Zs = self.calc_grid_impedence(pu=True, fl=new_tests["Grid_MVA_v"], x2r=new_tests["Grid_X2R_v"])
        if 'Vgrid' in row:
            ic("Vgrid in row")
            row_sect = {'Init_Vpoc_pu_v': self.calc_vpoc_from_vslack(vslack=row["Vgrid"], ppoc=ppoc, qpoc=qpoc, Zs=Zs),
                        'Vslack_pu_v': row["Vgrid"]}
        elif 'Vpoc' in row:
            ic("Vpoc in row")
            row_sect = {'Init_Vpoc_pu_v': row["Vpoc"],
                        'Vslack_pu_v': self.calc_vslack_from_vpoc(vpoc=row["Vpoc"], ppoc=ppoc, qpoc=qpoc, Zs=Zs)}
        else:
            row_sect = {'Init_Vpoc_pu_v': self.v_set,
                        'Vslack_pu_v': self.calc_vslack_from_vpoc(vpoc=self.v_set, ppoc=ppoc, qpoc=qpoc, Zs=Zs)}
        ic(row_sect)
        return [row_sect]
    
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
        return [row_sect]
    
    def add_fault_specs(self, row: pd.DataFrame, new_tests: list):
        ic("add_fault_specs")
        row_sect = dict()
        row_sect_list = []
        if 'Apply fault (s)' in row:
            ic("Apply fault")
            if 'Fault impedance (pu)' in row:
                ic("Fault impedance")
                if 'Fault duration (s)' in row:
                    ic("Fault duration")
                    if not type(row["Fault duration (s)"]) == str:
                        ic("is numeric")
                        if 'Fault type' in row:
                            ic("fault type")
                            x2r = new_tests["Grid_X2R_v"]
                            grid_impedance = self.calc_grid_impedence(pu=True, fl=new_tests["Grid_MVA_v"], x2r=x2r)
                            # Note this is not flexible: Only option for multiple fault impedances is "Zf=Rf = 1, 5 and 10 Ohm" and it must look like this exactly
                            ic(row["Fault impedance (pu)"])
                            
                            if row["Fault impedance (pu)"] == "Zf=Rf = 1, 5 and 10 Ohm":
                                fault_impedance_strs = ["Zf=Rf=1", "Zf=Rf=5", "Zf=Rf=10"]
                            else:
                                fault_impedance_strs = [row["Fault impedance (pu)"]]
                            ic(fault_impedance_strs)
                            for fault_impedance_str in fault_impedance_strs:
                                (fault_scr, fault_x2r, fault_impedance) = self.calc_fault_impedance(zf_str=fault_impedance_str, grid_impedance=grid_impedance, grid_x2r=x2r)
                                ic(fault_impedance_str)
                                ic(fault_scr)
                                Udip = self.calc_fault_voltage(Zs=grid_impedance, Zf=fault_impedance)
                                row_sect = {'Init_Fault_MVA': fault_scr*self.p_nom,
                                            'Init_Fault_X_on_R': fault_x2r,
                                            'Fault_Type_v': self.get_fault_code(row["Fault type"]),
                                            'Fault_Depth': Udip,
                                            'Fault_Time': row["Apply fault (s)"],
                                            'Fault_Duration': row["Fault duration (s)"],
                                            'Post_Fault_Duration': row["End Run (s)"] - row["Apply fault (s)"] - row["Fault duration (s)"],
                                            'Enable_Fault_Studies_v': 1,
                                            'Fault_Strategy_v': 0,
                                            'Fault_X2R_v': fault_x2r,
                                            'Fault_Timing_Signal_t': f'[0, {row["Apply fault (s)"]}, {row["Apply fault (s)"] + row["Fault duration (s)"]}]',
                                            'Fault_Timing_Signal_v': f'[0, 1, 0]',
                                            'Ures_v': Udip}
                        
                                row_sect_list.append(row_sect)
                        elif row["Fault impedance (pu)"][0:5] == "Yf = jXc":
                            
                        else:
                            raise CalcSheetError("fault type")
                    elif row["Fault duration (s)"] in ["S1", "S2", "S3", "S4", "S5"]:
                        row_sect = {'Enable_Fault_Studies_v': 1,
                                    'Fault_Strategy_v': 0,}
                        row_sect.update(self.mfrt_seq(row, new_tests))
                        row_sect_list.append(row_sect)
                    else:
                        raise CalcSheetError("fault duration")
                else:
                    raise CalcSheetError("fault duration")
            else:
                raise CalcSheetError("fault impedance")
        ic(row_sect_list)
        return row_sect_list
        
##################### MFRT GENERATOR #######################

    def mfrt_seq(self, row: pd.DataFrame, new_tests: list):
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
        row_sect.update({'Init_Fault_MVA': new_tests["Grid_SCR_v"]*self.p_nom,
                        'Init_Fault_X_on_R': new_tests["Grid_X2R_v"],
                        'Fault_X2R_v': new_tests["Grid_X2R_v"]})
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

#TODO    
    # def generate_profile(self, profile_name: str):
    #     figure_references = pd.read_excel()
    
    def calc_ramp(self, row: pd.DataFrame):
        ic('calc_ramp')
        d_freq = abs(row["Freq target"] - self.f_nom)
        d_t = d_freq/row["Freq ramp Hz/s"]
        t_end_ramp = min(row["Apply step (s)"] + d_t, row["End Run (s)"])
        return t_end_ramp
        
        
##################### FAULT CALCULATIONS #######################
# This is what we can read: "{number}", "Zf=0", "Zf=Zs", "Zf={number}xZs[Udip=~{number}pu]" (ignore [] part), "Zf=Rf=1, 5 and 10 Ohm", "Yf=jXc(U_Ov=1.15pu)" (calc Yf from U_Ov)
    def calc_fault_impedance(self, zf_str: str, grid_impedance: tuple, grid_x2r: float, ):
        ic('calc_fault_impedance')
        (_, _, grid_impedance) = grid_impedance
        # can't set fault impedance to 0 so set it to a very low number.
        min_fault_impedance = 0.000000001
        # can't set value to inf so set to very high value
        max_fault_x2r = 99999

        fault_x2r = grid_x2r
        ic(zf_str)
        if zf_str.isnumeric():
            if zf_str == 0:
                fault_impedance = min_fault_impedance
            else:
                fault_impedance = zf_str
        else:
            str_list = zf_str.split("=")
            ic(str_list)
            # If there is a "[comment]" part, remove this.
            if "[" in str_list[1]:
                str_list[1] = str_list[1].split("[")[0]
                str_list = str_list[0:2]
                ic(str_list)
            if len(str_list) >= 2 and str_list[0] == "Zf":
                if str_list[1] == "0":
                    # "Zf=0", str_list = [Zf, 0]
                    fault_impedance = min_fault_impedance
                elif str_list[1] == "Zs":
                    # "Zf=Zs", str_list = [Zf, Zs]
                    fault_impedance = grid_impedance
                elif str_list[1] == "Rf" and str_list[2].isnumeric():
                    # Zf=Rf=1
                    fault_impedance = float(str_list[2])
                    fault_x2r = 0
                else:
                    str_list = str_list[1].split("x")
                    if len(str_list) == 2 and str_list[1] == "Zs" and str_list[0].isnumeric():
                        # "Zf=2xZs" str_list = [2, Zs]
                        multiplier = int(str_list[0])
                        fault_impedance = grid_impedance * multiplier
                    else:
                        raise CalcSheetError("fault impedance")
            elif len(str_list) == 3 and str_list[0] == "Yf" and str_list[1] == "jXc(U_Ov":
                # Yf=jXc(U_Ov=1.15pu), str_list = ["Yf", "jXc(U_Ov", "1.15pu)"]
                str_list = str_list[2].split("pu")
                # str_list = [1.15, )]
                if len(str_list) == 2 and str_list[1] == ")":
                    u_ov = float(str_list[0])
                    fault_impedance = grid_impedance * u_ov/(1 + u_ov)
                    fault_x2r = max_fault_x2r
                else:
                    raise CalcSheetError("fault impedance")
            else:
                raise CalcSheetError("fault impedance")
        ic(self.v_nom)
        ic(self.p_nom)
        ic(fault_impedance)
        fault_scr = float(self.v_nom)**2/(fault_impedance*float(self.p_nom))
        # calculate a new fault impedence for every X/R ratio and SCR and return a list of tuples
        return (fault_scr, fault_x2r, fault_impedance)
    
    # def calc_fault_multiplier(self, string: str):
    #     multiplier = 0
    #     if string == "Zf=0":
    #         multiplier = 0
    #     elif string == "Zf=Zs":
    #         multiplier = 1
    #     return multiplier
    
    def calc_fault_voltage(self, Zf: tuple, Zs: tuple):
        ic('calc_fault_voltage')
        (_, _, zs) = Zs
        Udip= Zf/(zs + Zf)
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
        
##################### GENERAL CALCULATIONS #######################
    def read_scr_and_x2r(self, row: pd.DataFrame):
        ic("read_scr_and_x2r")
        if type(row["SCR"]) == str and row["SCR"][0] == "[" and row["SCR"][-1] == "]":
            scr_strs = json.loads(row["SCR"])
        else:
            scr_strs = [row["SCR"]]
        if type(row["X/R"]) == str and row["X/R"][0] == "[" and row["X/R"][0] == "]":
            x2r_strs = json.loads(row["X/R"])
        else:
            x2r_strs = [row["X/R"]]
        scr_and_x2r = []
        poc_scr = json.loads(self.system_inf.loc["POC SCR"]["Var Val"])
        poc_x2r = json.loads(self.system_inf.loc["POC XR ratio"]["Var Val"])
        for scr_str in scr_strs:
            for x2r_str in x2r_strs:
                if type(scr_str) == float or type(scr_str) == int:
                    if type(x2r_str) == float or type(x2r_str) == int:
                        scr_and_x2r.append((scr_str, x2r_str))
                    elif x2r_str == "POC":
                        scr_and_x2r.append((scr_str, poc_x2r[0]))
                        scr_and_x2r.append((scr_str, poc_x2r[1]))
                    else:
                        raise CalcSheetError("X2R")
                elif scr_str == "POC":
                    if type(x2r_str) == float or type(x2r_str) == int:
                        scr_and_x2r.append((poc_scr[0], x2r_str))
                        scr_and_x2r.append((poc_scr[1], x2r_str))
                    elif x2r_str == "POC":
                        scr_and_x2r.append((poc_scr[0], poc_x2r[0]))
                        scr_and_x2r.append((poc_scr[1], poc_x2r[1]))
                    else:
                        raise CalcSheetError("X2R")
                else:
                    raise CalcSheetError("SCR")
        return scr_and_x2r
    
    def calc_fault_level(self, scr: float):
        ic("calc_fault_level")
        fault_level = scr * self.p_nom
        ic(fault_level)
        return fault_level
    
    def calc_grid_impedence(self, pu: bool, fl: float, x2r: float):
        ic("calc_grid_impedence")
        z = self.v_nom**2 / fl
        r = z / math.sqrt(1+x2r**2)
        x = r * x2r
        if pu == True:
            r = r/self.z_base
            x = x/self.z_base
            z = z/self.z_base
        l = x / (2 * self.f_nom * math.pi)
        return (r, x, z)
        
    def calc_vslack_from_vpoc(self, vpoc: float, ppoc: float, qpoc: float, Zs: tuple):
        ic("calc_vslack_from_vpoc")
        spoc = math.sqrt(ppoc**2 + qpoc**2)
        (rgrid, xgrid, zgrid) = Zs
        rgrid_pu = rgrid/ self.z_base
        xgrid_pu = xgrid/ self.z_base
        zgrid_pu = zgrid/ self.z_base
        vslack = math.sqrt(vpoc**2 + spoc**2/vpoc**2*zgrid_pu**2 - 2*(ppoc*rgrid_pu + qpoc*xgrid_pu))
        return vslack
        
    def calc_vpoc_from_vslack(self, vslack: float, ppoc: float, qpoc: float, Zs: tuple):
        ic("calc_vpoc_from_vslack")
        spoc = math.sqrt(ppoc**2 + qpoc**2)
        (rgrid, xgrid, zgrid) = Zs
        rgrid_pu = rgrid/ self.z_base
        xgrid_pu = xgrid/ self.z_base
        zgrid_pu = zgrid/ self.z_base
        vgrid = math.sqrt((vslack**2 + 2*(ppoc*rgrid_pu + qpoc*xgrid_pu) + math.sqrt((vslack**2 + 2*(ppoc*rgrid_pu + qpoc*xgrid_pu))**2 - 4*spoc**2*zgrid_pu**2))/2)
        ic(vgrid)
        return vgrid
    
    