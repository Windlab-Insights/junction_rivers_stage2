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
        self.figure_references = pd.read_excel(calc_sheet_path, sheet_name="Figure References", index_col="Figure", header=0)
        
        # Record project information
        self.v_nom = self.system_inf.loc["POC base Voltage"]["Var Val"]
        self.f_nom = self.system_inf.loc["Nominal Frequency"]["Var Val"]
        self.p_nom = self.system_inf.loc["Nominal Power"]["Var Val"]
        self.q_nom = self.p_nom * 0.395
        self.s_nom = math.sqrt(self.p_nom**2 + self.q_nom**2)
        self.z_base = self.v_nom**2/self.p_nom
        self.q_set = 0
        self.v_set = self.system_inf.loc["V_POC"]["Var Val"]
        self.p_max_bess = 70
        self.qv_droop = self.system_inf.loc["Q-V Droop"]["Var Val"]

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
        # if there are no specs to add, then just return the test line as is
        if len(specs) == 0:
            return [in_row]
        # create a new list of tests
        out_rows = []
        # for multiple specs, we want to add new tests for each one
        for index, spec in enumerate(specs):
            # create a copy of the existing test to update without changing the original one
            ic(in_row)
            temp = deepcopy(in_row)
            ic(temp)
            ic(type(temp))
            temp.update({'File_Name': f'{in_row["File_Name"]}{index}'})
            temp.update(spec)
            out_rows.append(temp) 
        return out_rows          
    
    def add_new_row(self, row: pd.DataFrame, category: str):
        ic('add_new_row')
        # Initialise a dict with category information and a list of dicts with just the empty dict.
        first_test = dict()
        first_test.update({'File_Name': f'{category}_test_{row["Test"]}'.replace("-","_")})
        first_test.update({'DIR': category})
        first_test.update({'En_SMIB_init_v': 0})
        first_test.update({'Infinite_Grid_v': 0})
        new_tests = [first_test]
        
        file_infos = self.add_file_info(row)
        new_tests = self.update_rows(file_infos, first_test)
        ic(new_tests)
        
        out_rows = []
        for new_test in new_tests:
            q_specs = self.add_q_specs(row, new_test)
            out_rows.extend(self.update_rows(q_specs, new_test))
        new_tests = out_rows[:]
        ic(new_tests)
        
        out_rows = []
        for new_test in new_tests:
            p_ref_specs = self.add_p_ref_specs(row, self.BESS_PZERO)
            out_rows.extend(self.update_rows(p_ref_specs, new_test))
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
            vref_specs = self.add_vref_specs(row, new_test)
            out_rows.extend(self.update_rows(vref_specs, new_test))
        new_tests = out_rows[:]
        ic(new_tests)
        
        out_rows = []
        for new_test in new_tests:
            qref_specs = self.add_qref_specs(row, new_test)
            out_rows.extend(self.update_rows(qref_specs, new_test))
        new_tests = out_rows[:]
        ic(new_tests)
        
        out_rows = []
        for new_test in new_tests:
            osc_specs = self.add_osc_specs(row, new_test)
            out_rows.extend(self.update_rows(osc_specs, new_test))
        new_tests = out_rows[:]
        ic(new_tests)
        
        out_rows = []
        for new_test in new_tests:
            fault_specs = self.add_fault_specs(row, new_test)
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
                'Grid_SCR': scr,
                'Grid_X2R_v': x2r,
                'Grid_MVA_v': self.calc_fault_level(scr)
            })
            row_sect_list.append(temp_sect)
        ic(row_sect_list)
        return row_sect_list
    
    def add_p_ref_specs(self, row: pd.DataFrame, wf_state: int):
        ic('add_p_ref_specs')
        row_sect = dict()
        if 'Active Power (pu)' in row:
            if 'Time Steps (s)' in row and 'Pref_deltas (pu)' in row:
                pref_profile = self.Profile()
                pref_profile.read_profile(v_data=row["Pref_deltas (pu)"], t_data=row["Time Steps (s)"])
                pref_v = ([(row["Active Power"] + delta)*self.p_nom for delta in pref_profile.deltas])
                if wf_state == self.WTG_PZERO:
                    row_sect = {'Pref_Wind_MW_v': 0,
                                'Pref_BESS_MW_v': pref_v,
                                'Pref_BESS_MW_t': pref_profile.time_steps}
                elif wf_state == self.BESS_PMAX:
                    row_sect = {'Pref_Wind_MW_v': pref_v,
                                'Pref_BESS_MW_v': self.p_max_bess,
                                'Pref_Wind_MW_t': pref_profile.time_steps}
                elif wf_state == self.BESS_PZERO:
                    row_sect = {'Pref_Wind_MW_v': pref_v,
                                'Pref_BESS_MW_v': 0,
                                'Pref_Wind_MW_t': pref_profile.time_steps}
                elif wf_state == self.BESS_PMIN:
                    row_sect = {'Pref_Wind_MW_v': pref_v,
                                'Pref_BESS_MW_v': -self.p_max_bess,
                                'Pref_Wind_MW_t': pref_profile.time_steps}
                else:
                    raise CalcSheetError("wf_state")
            else:
                pref_v = row["Active Power (pu)"]*self.p_nom
                if wf_state == self.WTG_PZERO:
                    row_sect = {'Pref_Wind_MW_v': 0,
                                'Pref_BESS_MW_v': pref_v}
                elif wf_state == self.BESS_PMAX:
                    row_sect = {'Pref_Wind_MW_v': pref_v,
                                'Pref_BESS_MW_v': self.p_max_bess}
                elif wf_state == self.BESS_PZERO:
                    row_sect = {'Pref_Wind_MW_v': pref_v,
                                'Pref_BESS_MW_v': 0}
                elif wf_state == self.BESS_PMIN:
                    row_sect = {'Pref_Wind_MW_v': pref_v,
                                'Pref_BESS_MW_v': -self.p_max_bess}
                else:
                    raise CalcSheetError("wf_state")
        else:
            row_sect = {'Pref_Wind_MW_v': 0,
                        'Pref_BESS_MW_v': 0}
        ic(row_sect)
        return [row_sect]
    
    def add_q_specs(self, row: pd.DataFrame, new_tests: list):
        ic("add_q_specs")
        row_sect = dict()
        # Reactive power is specified
        if 'Reactive Power (pu)' in row:
            row_sect = {'Init_Qpoc_pu_v': row["Reactive Power (pu)"] * self.s_nom}
        else:
            row_sect = {'Init_Qpoc_pu_v': 0}
        return [row_sect]
            
    def add_v_specs(self, row: pd.DataFrame, new_tests: list):
        ic('add_v_specs')
        row_sect_list = []
        row_sect = dict()
        # Use data which has already been added to this test to calculate the values to be added
        # If the active power reference is a list, then use the first one as a reference
        if type(new_tests["Pref_Wind_MW_v"]) == list:
            pref_wind = new_tests["Pref_Wind_MW_v"][0]
        else:
            pref_wind = new_tests["Pref_Wind_MW_v"]
        if type(new_tests["Pref_BESS_MW_v"]) == list:
            pref_bess = new_tests["Pref_BESS_MW_v"][0]
        else:
            pref_bess = new_tests["Pref_BESS_MW_v"]
        ppoc = pref_wind + pref_bess
        qpoc = new_tests["Init_Qpoc_pu_v"]
        Zs = self.calc_grid_impedence(pu=True, fl=new_tests["Grid_MVA_v"], x2r=new_tests["Grid_X2R_v"])
        
        # Vgrid/Vslack is specified, calculate V POC
        if 'Vgrid' in row:
            row_sect = {'Init_Vpoc_pu_v': self.calc_vpoc_from_vslack(vslack=row["Vgrid"], ppoc=ppoc, qpoc=qpoc, Zs=Zs),
                        'Vslack_pu_v': row["Vgrid"]}
            row_sect_list = [row_sect]
        # Vpoc is specified, calculate Vgrid/Vslack
        elif 'Vpoc' in row:
            row_sect = {'Init_Vpoc_pu_v': row["Vpoc"],
                        'Vslack_pu_v': self.calc_vslack_from_vpoc(vpoc=row["Vpoc"], ppoc=ppoc, qpoc=qpoc, Zs=Zs)}
            row_sect_list = [row_sect]
        # There is a Vgrid profile in DMAT tests 149 - 166
        elif 'Event' in row and row["Event"] == "Vgrid":
            vgrid_profile = self.Profile()
            vgrid_profile.read_profile(row)
            # vslack starting point is calculated form Vpoc
            vslack_starting_point = self.calc_vslack_from_vpoc(vpoc=self.v_set, ppoc=ppoc, qpoc=qpoc, Zs=Zs)
            row_sect = {'Init_Vpoc_pu_v': self.v_set,
                        'Vslack_pu_v': [delta + vslack_starting_point for delta in vgrid_profile.deltas],
                        'Vslack_pu_t': vgrid_profile.time_steps}
            row_sect_list = [row_sect]
        elif 'Volt_delta (pu)' in row and 'Time Steps (s)' in row and 'Ramp/Step' in row:
            vgrid_profile = self.Profile()
            vslack_starting_point = self.calc_vslack_from_vpoc(vpoc=self.v_set, ppoc=ppoc, qpoc=qpoc, Zs=Zs)
            if vgrid_profile.check_valid_input(row["Volt_delta (pu)"]):
                vgrid_profile.read_profile(vdata=row["Volt_delta (pu)"], t_data=row["Time Steps (s)"], r_data=row["Ramp/Step"])
                row_sect = {'Init_Vpoc_pu_v': self.v_set,
                            'Vslack_pu_v': [delta + vslack_starting_point for delta in vgrid_profile.deltas],
                            'Vslack_pu_t': vgrid_profile.time_steps,
                            'Vslack_pu_r': vgrid_profile.ramps}
                row_sect_list = [row_sect]
            else:
                figures = row["Volt_delta (pu)"].split(",")
                for figure in figures:
                    if vgrid_profile.check_valid_input(row["Volt_delta (pu)"]):
                        vgrid_profile.read_profile(fig=figure)
                        row_sect = {'Init_Vpoc_pu_v': self.v_set,
                                    'Vslack_pu_v': [delta + vslack_starting_point for delta in vgrid_profile.deltas],
                                    'Vslack_pu_t': vgrid_profile.time_steps,
                                    'Vslack_pu_r': vgrid_profile.ramps}
                        row_sect_list.append(row_sect)
                    else:
                        raise CalcSheetError("volt delta")
        else:
            row_sect = {'Init_Vpoc_pu_v': self.v_set,
                        'Vslack_pu_v': self.calc_vslack_from_vpoc(vpoc=self.v_set, ppoc=ppoc, qpoc=qpoc, Zs=Zs)}
            row_sect_list = [row_sect]
        ic(row_sect)
        return row_sect_list
    
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
    
    def add_qref_specs(self, row: pd.DataFrame, new_tests: list):
        ic("add_qref_specs")
        row_sect = dict()
        if 'Event' in row and row["Event"] == "Qref":
            qref_profile = self.Profile()
            qref_profile.read_profile(v_data="Delta (pu)", t_data="Time Steps (s)")
            q_starting_point = new_tests["Init_Qpoc_pu_v"]
            row_sect = {'Qref_MVAr_v': [(delta + q_starting_point) * self.q_nom for delta in qref_profile.deltas],
                        'Qref_MVAr_t': qref_profile.time_steps}
        return [row_sect]
    
    def add_vref_specs(self, row: pd.DataFrame, new_tests: list):
        ic("add_vref_specs")
        # Get the Q ref value in case it is needed
        if 'Qref_MVAr_v' in new_tests:
            qref = new_tests["Qref_MVAr_v"][0]
        else:
            qref = self.q_set
        
        row_sect = dict()
        if 'Event' in row and row["Event"] == "Vref":
            vref_profile = self.Profile()
            vref_profile.read_profile(v_data="Delta (pu)", t_data="Time Steps (s)")
            v_starting_point = new_tests["Init_Vpoc_pu_v"]
            row_sect = {'Vref_pu_v': [delta + v_starting_point for delta in vref_profile.deltas],
                        'Vref_pu_t': vref_profile.time_steps}
        # We might want to calculate Vref based off the given Qpoc and Vpoc values
        elif 'Init_Vpoc_pu_v' in new_tests and 'Init_Qpoc_MVAr_v' in new_tests:
            vpoc = new_tests["Init_Vpoc_pu_v"]
            qpoc = new_tests["Init_Qpoc_MVAr_v"]
            row_sect = {'Vref_pu_v': self.calc_vref_from_qpoc_and_vpoc(qpoc=qpoc,vpoc=vpoc,qref=qref)}
        else:
            row_sect = {'Vref_pu_v': self.v_set}
        # # Update the qpoc value based off vref and vpoc    
        # if 'Init_Vpoc_pu_v' in new_tests and not 'Init_Qpoc_MVAr_v' in new_tests:
        #     vpoc = new_tests["Init_Vpoc_pu_v"]
        #     vref = row_sect["Vref_pu_v"][0]
        #     row_sect.update({'Init_Qpoc_MVAr_v': self.calc_qpoc_from_vref_and_vpoc(vref=vref,vpoc=vpoc,qref=qref)})
        # # Update the vpoc value based off vref and qpoc
        # elif not 'Init_Vpoc_pu_v' in new_tests and 'Init_Qpoc_MVAr_v' in new_tests:
        #     qpoc = new_tests["Init_Qpoc_MVAr_v"]
        #     vref = row_sect["Vref_pu_v"][0]
        #     row_sect.update({'Init_Vpoc_pu_v': self.calc_vpoc_from_qpoc_and_vref(qpoc=qpoc,vref=vref,qref=qref)})
            
        return [row_sect]
    
    def add_osc_specs(self, row: pd.DataFrame, new_tests: list):
        row_sect_list = []
        row_sect = dict()
        if 'Time Steps (s)' in row:
            if 'Osc_freq' in row:
                if 'Magnitude' in row:
                    if 'Step (Hz)' in row:
                        if 'Osc_Phase_deg' in row:
                            # read the time steps
                            try:
                                time_steps = json.loads(row["Time Steps (s)"])
                            except:
                                raise CalcSheetError("time steps")
                            timing_sig = [0, time_steps[0], time_steps[1], row["End Run (s)"]]
                            # split the freq steps and osc freqs into different test sets
                            freq_steps_strs = str(row["Step (Hz)"]).split(",")
                            osc_freq_strs = str(row["Osc_freq"]).split(",")
                            # check that the number of test groups is the same and raise an error if not
                            if not len(osc_freq_strs) == len(freq_steps_strs):
                                raise CalcSheetError("time steps and freq steps")
                            # iterrate through each test group
                            for index, osc_freq_str in enumerate(osc_freq_strs):
                                osc_freqs = osc_freq_str.split(":")
                                step = float(freq_steps_strs[index])
                                # make a list of all the frequencies which should be added
                                osc_freq_range = []
                                temp = float(osc_freqs[0])
                                while temp < float(osc_freqs[1]):
                                    osc_freq_range.append(temp)
                                    temp = temp + step
                                # itterate through each frequency level and add a test
                                for osc_freq in osc_freq_range:
                                    row_sect = {"Vslack_osc_amplitude_v": [0, row["Magnitude"], 0, 0],
                                                "Vslack_osc_amplitude_t": timing_sig,
                                                "Vslack_osc_Hz_v": [0, osc_freq, 0, 0],
                                                "Vslack_osc_Hz_t": timing_sig,
                                                "Vslack_osc_phase_deg_v": [0, row["Osc_Phase_deg"], 0, 0],
                                                "Vslack_osc_phase_deg_t": timing_sig}
                                    row_sect_list.append(row_sect)
                        else:
                            raise CalcSheetError("osc phase angle")
                    else:
                        raise CalcSheetError("step")
                else:
                    raise CalcSheetError("magnitude")
            elif 'Magnitude' in row or 'Osc_Phase_deg' in row:
                raise CalcSheetError('Osc frequency')
        return row_sect_list
    
    def is_number(self, string:str):
        try:
            float(string)
            return True
        except ValueError:
            return False
    
    def add_fault_specs(self, row: pd.DataFrame, new_tests: list):
        ic("add_fault_specs")
        row_sect = dict()
        row_sect_list = []
        if 'Apply fault (s)' in row:
            if 'Fault impedance (pu)' in row:
                if 'Fault duration (s)' in row:
                    fault_duration = str(row["Fault duration (s)"])
                    if self.is_number(fault_duration):
                        if 'Fault type' in row:
                            x2r = new_tests["Grid_X2R_v"]
                            ic(row["Fault impedance (pu)"])
                            # TODO: not sure what Enable fault studies means
                            # Fault strategy=1 means that the "Zs2Zf method" is used which means that the impedance is calculated as a ratio of the grid impedance.
                            # The alternative is fault startegy =0 which means the fault impedance is calculated based off Udip in PSCAD.
                            row_sect = {'Enable_Fault_Studies_v': 1,
                                        'Fault_Strategy_v': 1,
                                        'Fault_X2R_v': x2r,
                                        'Fault_Type_v': self.get_fault_code(row["Fault type"]),
                                        'Fault_Timing_Signal_t': f'[0, {row["Apply fault (s)"]}, {row["Apply fault (s)"] + row["Fault duration (s)"]}]',
                                        'Fault_Timing_Signal_v': f'[0, 1, 0]'}
                            # Resistive faults
                            if row["Fault impedance (pu)"] == "Zf=Rf = 1, 5 and 10 Ohm":
                            # Note this is not flexible: Only option for multiple fault impedances is "Zf=Rf = 1, 5 and 10 Ohm" and it must look like this exactly
                                R_ohms = [1, 5, 10]
                                for R_ohm in R_ohms:
                                    row_sect_copy = deepcopy(row_sect)
                                    row_sect_copy.update({'Zf2Zs_v': 0, 
                                                     'Rf_Offset_ohms_v': R_ohm,
                                                     'Xf_Offset_ohms_v': 0})
                                    row_sect_list.append(row_sect_copy)
                            # Fault impedance specified like "Zf=0", "Zf=Zs", or "Zf=2xZf"
                            else:
                                row_sect.update({'Zf2Zs_v': self.read_fault_multiplier(row["Fault impedance (pu)"])})
                                row_sect_list.append(row_sect)
                        # Temporary over voltage specfied as "Yf = jXc (U_Ov = 1.15pu)"
                        elif "Yf = jXc" in row["Fault impedance (pu)"] and "U_Ov =" in row["Fault impedance (pu)"]:
                            # 'Yf = jXc (U_Ov = 1.15pu)' --> ['Yf', 'jXc (U_Ov', '1.15pu') ]
                            fault_voltage = row["Fault impedance (pu)"].split("=")
                            # '1.15pu)' --> ['1.15', ')']
                            fault_voltage = fault_voltage[2].split("pu")[0]
                            grid_impedance = self.calc_grid_impedence(pu=1, fl=new_tests["Grid_MVA_v"], x2r=new_tests["Grid_X2R_v"])
                            fault_impedance = self.calc_fault_impedance(fault_voltage=float(fault_voltage), fault_distance=1, grid_impedance=grid_impedance)
                            row_sect = {'TOV_Shunt_uF_v': 1/(fault_impedance * self.f_nom * 2 * math.pi) * 10**6,
                                        'TOV_Timing_Signal_t': f'[0, {row["Apply fault (s)"]}, {row["Apply fault (s)"] + row["Fault duration (s)"]}]',
                                        'TOV_Timing_Signal_v': f'[0, 1, 0]'}
                            row_sect_list = [row_sect]
                    elif fault_duration in ["S1", "S2", "S3", "S4", "S5"]:
                        ic(row["Fault duration (s)"])
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
        random.seed(row["Test"])
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
        row_sect.update({'Fault_X2R_v': new_tests["Grid_X2R_v"]})
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
        
        row_sect = {'Fault_Timing_Signal_t': fault_timing_signal_t,
                    'Fault_Timing_Signal_v': fault_timing_signal_v,
                    'Zf2Zs_t': fault_times,
                    'Zf2Zs_v': fault_zs_multiplier,
                    'Fault_Type_t': fault_times,
                    'Fault_Type_v': fault_types}    
        return row_sect
    
    # To simplify access to a profile, there is a profile subclass which holds the information for a signal including _v, _t, and _r information
    class Profile:
        def __init__(self):
            self.deltas = None
            self.time_steps = None
            self.ramps = None
            self.profile_type = None
    
        def check_valid_input(self, string):
            try:
                _ = json.loads(string)
                return True
            except:
                if string in self.figure_references["Figure"]:
                    return True
                else:
                    return False
        
        def get_default_profile_type(self, figure_name: str):
            figure_reference = self.figure_references.loc[figure_name]
            self.profile_type =  figure_reference["Type"]
            
        def read_profile(self, fig: str):
            # default profile is refered to as "fig"
            try:
                figure_reference = self.figure_references.loc[fig]
                self.deltas = json.loads(figure_reference["Deltas"])
                self.time_steps = json.loads(figure_reference["Time_steps"])
                self.get_default_profile_type(fig)
            except:
                raise CalcSheetError("figure name")
        
        def read_profile(self, v_data: str, t_data: str, r_data=None):
            # manual profile is refered to as [1,2,3]
            try:
                self.deltas = json.loads(v_data)
                self.time_steps = json.loads(t_data)
                if not r_data == None:
                    try:
                        self.ramps = json.loads(r_data)
                    except:
                        raise CalcSheetError("ramps")
            except:
                try:
                    self.read_profile(v_data)
                except:
                    raise CalcSheetError("time steps or delta")
                
    
    def calc_ramp(self, row: pd.DataFrame):
        ic('calc_ramp')
        d_freq = abs(row["Freq target"] - self.f_nom)
        d_t = d_freq/row["Freq ramp Hz/s"]
        t_end_ramp = min(row["Apply step (s)"] + d_t, row["End Run (s)"])
        return t_end_ramp
    
    def generate_oscillation(self, vslack: float, ):
        ic("generate_oscillation")
        
        
        
##################### FAULT CALCULATIONS #######################
    def calc_fault_impedance(self, fault_voltage: float, fault_distance: float, grid_impedance: tuple):
        ic('calc_fault_impedance')
        (_, _, grid_impedance) = grid_impedance
        ic(grid_impedance)
        ic(fault_voltage)
        grid_impedance = grid_impedance*self.z_base
        fault_impedance = fault_distance*grid_impedance*fault_voltage/abs(1 - fault_voltage)
            
        return fault_impedance
    
    def read_fault_multiplier(self, zf_str: str):
        ic('read_fault_multiplier')
        # can't set fault impedance to 0 so set it to a very low number.
        min_fault_impedance = 0.000000001
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
                multiplier = min_fault_impedance
            elif str_list[1] == "Zs":
                # "Zf=Zs", str_list = [Zf, Zs]
                multiplier = 1
            else:
                # "Zf=2xZs", str_list = [Zf, 2xZs]
                str_list = str_list[1].split("x")
                if len(str_list) == 2 and str_list[1] == "Zs" and str_list[0].isnumeric():
                    # "2xZs" str_list = [2, Zs]
                    multiplier = int(str_list[0])
                else:
                    raise CalcSheetError("fault impedance")
        else:
            raise CalcSheetError("fault impedance")
        return multiplier
    
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
        scr_strs = str(row["SCR"]).split(",")
        x2r_strs = str(row["X/R"]).split(",")
        # if type(row["SCR"]) == str and row["SCR"][0] == "[" and row["SCR"][-1] == "]":
        #     scr_strs = json.loads(row["SCR"])
        # else:
        #     scr_strs = [row["SCR"]]
        # if type(row["X/R"]) == str and row["X/R"][0] == "[" and row["X/R"][-1] == "]":
        #     x2r_strs = json.loads(row["X/R"])
        # else:
        #     x2r_strs = [row["X/R"]]
        scr_and_x2r = []
        poc_scr = json.loads(self.system_inf.loc["POC SCR"]["Var Val"])
        poc_x2r = json.loads(self.system_inf.loc["POC XR ratio"]["Var Val"])
        for scr_str in scr_strs:
            for x2r_str in x2r_strs:
                if scr_str.isnumeric():
                    if x2r_str.isnumeric():
                        scr_and_x2r.append((float(scr_str), float(x2r_str)))
                    elif x2r_str == "POC":
                        scr_and_x2r.append((float(scr_str), float(poc_x2r[0])))
                        scr_and_x2r.append((float(scr_str), float(poc_x2r[1])))
                    else:
                        raise CalcSheetError("X2R")
                elif scr_str == "POC":
                    if x2r_str.isnumeric():
                        scr_and_x2r.append((float(poc_scr[0]), float(x2r_str)))
                        scr_and_x2r.append((float(poc_scr[1]), float(x2r_str)))
                    elif x2r_str == "POC":
                        scr_and_x2r.append((float(poc_scr[0]), float(poc_x2r[0])))
                        scr_and_x2r.append((float(poc_scr[1]), float(poc_x2r[1])))
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
    
    def calc_qpoc_from_vref_and_vpoc(self, vref: float, vpoc: float, qref: float):
        ic("calc_qpoc_from_vref_and_vpoc")
        qpoc = (vpoc - vref)/self.qv_droop + qref
        return qpoc
    
    def calc_vref_from_qpoc_and_vpoc(self, qpoc: float, vpoc: float, qref: float):
        ic("calc_vref_from_qpoc_and_vpoc")
        vref = vpoc - (qpoc - qref)*self.qv_droop
        return vref
    
    def calc_vpoc_from_qpoc_and_vref(self, qpoc: float, vref: float, qref: float):
        ic("calc_vpoc_from_qpoc_and_vref")
        vpoc = (qpoc - qref)*self.qv_droop + vref
        return vpoc