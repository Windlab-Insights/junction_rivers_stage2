import pandas as pd
import json
from typing import Dict, Optional
import math
from icecream import ic
import random
from copy import deepcopy
from enum import Enum
import os
from rengen.utils.time_utils import get_date_time_str

class CalcSheetError(Exception):
    def __init__(self):
        print(f"Calc Sheet Error")
    
    def __init__(self, input):
        print(f"Calc Sheet Error: Cannot interpret {input}")
        

class SpecGenerator():
    
    class WfState(Enum):
     BESS_DMAT_WTG_PZERO = 0
     WTG_DMAT_BESS_PMAX = 1
     WTG_DMAT_BESS_PZERO = 2
     WTG_DMAT_BESS_PMIN = 3
     BESS_PNEG_DMAT_WTG_PZERO = 4
    
    HEADERS = [
        'Grouping',
        'Key_Tests',
        'File_Name',
        'Post_Init_Duration_s',
        'Category',
        'DIR',
        ]
    
    def __init__(self, config_sheet_path, spec_output_dir, suite, run_all_tests):
        self.spec_generator(config_sheet_path, spec_output_dir, suite, run_all_tests)
        
    def spec_generator (self, config_sheet_path, spec_output_dir, suite, run_all_tests):
        date_time = get_date_time_str()
        # Read in spec sheet info from the Config sheet as data frames
        self.system_inf = pd.read_excel(config_sheet_path, sheet_name="System_inf", index_col="Var name", header=0)
        ic(self.system_inf)
        # Record project information
        self.v_nom = self.system_inf.loc["POC Base Voltage"]["Var Val"]
        self.f_nom = self.system_inf.loc["Nominal Frequency"]["Var Val"]
        self.p_nom = self.system_inf.loc["Nominal Active Power"]["Var Val"]
        self.q_nom = self.system_inf.loc["Nominal Reactive Power"]["Var Val"]
        self.s_nom = math.sqrt(self.p_nom**2 + self.q_nom**2)
        self.z_base = self.v_nom**2/self.p_nom
        self.q_set = 0
        self.v_set = self.system_inf.loc["Nominal Voltage"]["Var Val"]
        self.p_max_bess = self.system_inf.loc["BESS Capacity"]["Var Val"]
        self.qv_droop = self.system_inf.loc["Q-V Droop"]["Var Val"]
        self.software = self.system_inf.loc["Software"]["Var Val"]
        self.vgrid_change_method = self.system_inf.loc["Vgrid"]["Var Val"]
        self.seperate_spec_sheets = self.system_inf.loc["Seperate Spec Sheets"]["Var Val"]
        
        # if a suite is specified then, access the relevant excel document to access the tests
        if suite == "csr":
            calc_sheet_path = self.system_inf.loc["CSR address"]["Var Val"]
        elif suite == "dmat" or suite == "r_dmat":
            calc_sheet_path = self.system_inf.loc["DMAT address"]["Var Val"]
        else:
            calc_sheet_path = config_sheet_path
            
        # Read in info from config sheet
        checklist = pd.read_excel(calc_sheet_path, sheet_name="Checklist", index_col="Checklist", header=0)
        self.figure_references = pd.read_excel(calc_sheet_path, sheet_name="Figure References", index_col="Figure", header=0)
        
        # Initialise data frame for the output spec data
        spec_df = pd.DataFrame()
        
        # Iterrate through the categories as recorded in the checklist
        for index_cat, category in checklist.iterrows():
            if category.loc["Action"] == "Yes":
                ic(category)
                category_formated = index_cat.lower().replace(" ", "_")
                tests = pd.read_excel(calc_sheet_path, sheet_name=index_cat, header=0)
                
                # Iterate through the tests in the sheet
                for _, row in tests.iterrows():
                    run_test = run_all_tests or row["Action"] == "Yes"
                    # Determine if to run the test or not based on the suite selection
                    if suite == "r_dmat" and run_test:
                        ic(row["Suite"] == "Reduced_DMAT")
                        run_test = row["Suite"] == "Reduced_DMAT"
                    # Determine if to run the test based on the software which the test is designed for
                    if self.software == "PSCAD":
                        run_test = run_test and (row["Software"] == "PSCAD" or row["Software"] == "Both")
                    elif self.software == "PSSE":
                        run_test = run_test and (row["Software"] == "PSSE" or row["Software"] == "Both")
                    # add rows to the spec sheet to actually run the test
                    if run_test:
                        ic(row)
                        new_row = self.add_new_row(row, category=category_formated)
                        for row in new_row:
                            row.update({'File_Name': f'{row["File_Name"]}_{len(spec_df.index)+2:03d}'})
                            row.update({'Category': category_formated})
                            spec_df = spec_df.append(row, ignore_index=True)
                if self.seperate_spec_sheets:
                    spec_path = os.path.join(spec_output_dir,f"{date_time}_{category_formated}.csv")
                    spec_df.to_csv(spec_path, index=False)
                    ic(spec_path)
                    spec_df = pd.DataFrame()
        if not self.seperate_spec_sheets:
            spec_path = os.path.join(spec_output_dir,f"{date_time}_runner_spec.csv")
            spec_df.to_csv(spec_path, index=False)
            ic(spec_path)
        
################# ADD PARAMS FOR A TEST ################
    # for rows which have just been added to the new spec sheet, add new information and split into new rows if needed
    def update_rows(self, specs: list, in_row: dict):
        # if there are no specs to add, then just return the test line as is
        if len(specs) == 0:
            return [in_row]
        # create a new list of tests
        out_rows = []
        # for multiple specs, we want to add new tests for each one
        for spec in specs:
            # create a copy of the existing test to update without changing the original one
            temp = deepcopy(in_row)
            temp.update(spec)
            out_rows.append(temp)
        return out_rows
    
    # for each test in the config file to be added, add the test information one at a time.
    def add_new_row(self, row: pd.DataFrame, category: str):
        # Initialise a dict with category information and a list of dicts with just the empty dict.
        first_test = dict()
        # remove columns where the cell is blank
        row = row.dropna()
        first_test.update({'File_Name': row["Test Name"]})
        first_test.update({'DIR': category})
        first_test.update({'Test_Number': f'{row["Test Name"]}_{row["Test Number"]}'})
        first_test.update({'En_SMIB_init_v': 0})
        new_tests = [first_test]
        
        file_infos = self.add_file_info(row)
        new_tests = self.update_rows(file_infos, first_test)
        
        # for each spec sheet row which has been added (a new run), add the new information.
        out_rows = []
        for new_test in new_tests:
            q_specs = self.add_q_specs(row, new_test)
            out_rows.extend(self.update_rows(q_specs, new_test))
            new_tests = out_rows[:]
        
        out_rows = []
        for new_test in new_tests:
            p_ref_specs = self.add_p_ref_specs(row, new_test)
            out_rows.extend(self.update_rows(p_ref_specs, new_test))
            new_tests = out_rows[:]
        
        out_rows = []
        for new_test in new_tests:
            v_specs = self.add_v_specs(row, new_test)
            out_rows.extend(self.update_rows(v_specs, new_test))
        new_tests = out_rows[:]
        
        out_rows = []
        for new_test in new_tests:
            freq_specs = self.add_freq_specs(row)
            out_rows.extend(self.update_rows(freq_specs, new_test))
        new_tests = out_rows[:]
        
        out_rows = []
        for new_test in new_tests:
            vref_specs = self.add_vref_specs(row, new_test)
            out_rows.extend(self.update_rows(vref_specs, new_test))
        new_tests = out_rows[:]
        
        out_rows = []
        for new_test in new_tests:
            qref_specs = self.add_qref_specs(row, new_test)
            out_rows.extend(self.update_rows(qref_specs, new_test))
        new_tests = out_rows[:]
        
        out_rows = []
        for new_test in new_tests:
            osc_specs = self.add_osc_specs(row, new_test)
            out_rows.extend(self.update_rows(osc_specs, new_test))
        new_tests = out_rows[:]
        
        out_rows = []
        for new_test in new_tests:
            phase_specs = self.add_phase_specs(row, new_tests)
            out_rows.extend(self.update_rows(phase_specs, new_test))
        new_tests = out_rows[:]
        
        out_rows = []
        for new_test in new_tests:
            fault_specs = self.add_fault_specs(row, new_test)
            out_rows.extend(self.update_rows(fault_specs, new_test))
        new_tests = out_rows[:]
        return new_tests
    
    # add general test information including run time, scr, and x/r levels.
    def add_file_info(self, row: pd.DataFrame):
        row_sect_list = []
      
        row_sect = dict()
        row_sect =  {
                    'Grouping': '',
                    'Key_Tests': '',
                    'Category': '',
                    'Post_Init_Duration_s': row["End Run (s)"],
                }
        if str(row["SCR"]) == "inf" and str(row["X/R"]) == "inf":
            row_sect.update({'Infinite_Grid_v': 1})
            row_sect_list = [row_sect]
        else:
            row_sect.update({'Infinite_Grid_v': 0})
            scr_and_x2rs = self.read_scr_and_x2r(row["SCR"], row["X/R"])
            for scr_and_x2r in scr_and_x2rs:
                (scr, x2r) = scr_and_x2r
                temp_sect = row_sect.copy()
                try:
                    scr_vals= json.loads(scr) # If this works, then we have multiple values of SCR in a single test
                    if 'Apply Fault (s)' in row:
                        temp_sect.update({
                            'Grid_SCR': scr_vals,
                            'Grid_X2R_v': x2r,
                            'Grid_MVA_v': [self.calc_fault_level(scr_val) for scr_val in scr_vals],
                            'Grid_MVA_t': [0, row["Apply Fault (s)"] + row["Fault Duration (s)"], row["End Run (s)"]],
                        })
                    else:
                        raise CalcSheetError("SCRs")
                except TypeError:
                    temp_sect.update({
                        'Grid_SCR': scr,
                        'Grid_X2R_v': x2r,
                        'Grid_MVA_v': self.calc_fault_level(scr)
                    })
                finally:
                    row_sect_list.append(temp_sect)
        return row_sect_list
    
    # Add Pref_Wind_MW and Pref_BESS_MW info
    def add_p_ref_specs(self, row: pd.DataFrame, new_test: dict):
        # from the system inf sheet, get the wf states and make a list of them to run through
        if self.system_inf.loc["WF States"]["Var Val"] == "ALL":
            wf_states = list(self.WfState)
        else:
            wf_state_strs = self.system_inf.loc["WF States"]["Var Val"].split("; ")
            wf_states = []
            for wf_state_str in wf_state_strs:
                if wf_state_str == "WTG PZERO":
                    wf_states.append(self.WfState.BESS_DMAT_WTG_PZERO)
                elif wf_state_str == "BESS PMAX":
                    wf_states.append(self.WfState.WTG_DMAT_BESS_PMAX)
                elif wf_state_str == "BESS PZERO":
                    wf_states.append(self.WfState.WTG_DMAT_BESS_PZERO)
                elif wf_state_str == "BESS PMIN":
                    wf_states.append(self.WfState.WTG_DMAT_BESS_PMIN)
                elif wf_state_str == "BESS PNEG": ##### consider what we want this to be called.
                    wf_states.append(self.WfState.BESS_PNEG_DMAT_WTG_PZERO)
                else:
                    raise CalcSheetError("wf state")
        row_sect = dict()
        row_sect_list = []
        if 'WTG Pref (pu)' in row and 'BESS Pref (pu)' in row:
            # WTG pref information
            if 'WTG Pref Deltas (pu)' in row and 'Time Steps (s)' in row:
                wtg_pref_profile = self.Profile(self.figure_references)
                if 'WTG Pref Ramp (pu/s)' in row:
                    wtg_pref_profile.read_profile(v_data=row["WTG Pref Deltas (pu)"], t_data=row["Time Steps (s)"], r_data=row["WTG Pref Ramp (pu/s)"])
                    row_sect = {'Pref_Wind_MW_v': [(row["WTG Pref (pu)"] + delta)*self.p_nom for delta in wtg_pref_profile.deltas],
                                'Pref_Wind_MW_t': wtg_pref_profile.time_steps,
                                'Pref_Wind_MW_r': wtg_pref_profile.ramps}
                else:
                    wtg_pref_profile.read_profile(v_data=row["WTG Pref Deltas (pu)"], t_data=row["Time Steps (s)"])
                    row_sect = {'Pref_Wind_MW_v': [(row["WTG Pref (pu)"] + delta)*self.p_nom for delta in wtg_pref_profile.deltas],
                                'Pref_Wind_MW_t': wtg_pref_profile.time_steps}
            elif 'WTG Pref Deltas (pu)' in row:
                raise CalcSheetError("wtg pref time steps")
            else:
                row_sect = {'Pref_Wind_MW_v': float(row["WTG Pref (pu)"])*self.p_nom}
            # add BESS pref information
            if 'BESS Pref Deltas (pu)' in row and 'Time Steps (s)' in row:
                bess_pref_profile = self.Profile(self.figure_references)
                if 'BESS Pref Ramp (pu/s)' in row:
                    bess_pref_profile.read_profile(v_data=row["BESS Pref Deltas (pu)"], t_data=row["Time Steps (s)"], r_data=row["BESS Pref Ramp (pu/s)"])
                    row_sect.update({'Pref_BESS_MW_v': [(row["BESS Pref (pu)"] + delta)*self.p_max_bess for delta in bess_pref_profile.deltas],
                                    'Pref_BESS_MW_t': bess_pref_profile.time_steps,
                                    'Pref_BESS_MW_r': bess_pref_profile.ramps})
                else:
                    bess_pref_profile.read_profile(v_data=row["BESS Pref Deltas (pu)"], t_data=row["Time Steps (s)"])
                    row_sect.update({'Pref_BESS_MW_v': [(row["BESS Pref (pu)"] + delta)*self.p_max_bess for delta in bess_pref_profile.deltas],
                                    'Pref_BESS_MW_t': bess_pref_profile.time_steps})
            elif 'BESS Pref Deltas (pu)' in row:
                raise CalcSheetError("bess pref time steps")
            else:
                row_sect.update({'Pref_BESS_MW_v': float(row["BESS Pref (pu)"])*self.p_max_bess})
            # convert to list
            row_sect_list = [row_sect]
        elif 'Active Power (pu)' in row:
            for wf_state in wf_states:
                # scale the reference by the BESS pmax or WT pmax depending on the state
                if wf_state == self.WfState.BESS_DMAT_WTG_PZERO or wf_state == self.WfState.BESS_PNEG_DMAT_WTG_PZERO:
                    p_base = self.p_max_bess
                else:
                    p_base = self.p_nom
                # set a ceiling to ensure we don't go over the maximum P output of the wind farm
                if wf_state == self.WfState.WTG_DMAT_BESS_PMAX:
                    p_ceil = self.p_nom - self.p_max_bess
                else:
                    p_ceil = self.p_nom
                # add profile
                if 'Time Steps (s)' in row and 'Pref Deltas (pu)' in row:
                    pref_profile = self.Profile(self.figure_references)
                    pref_profile.read_profile(v_data=row["Pref Deltas (pu)"], t_data=row["Time Steps (s)"])
                    pref_v = ([min((row["Active Power (pu)"] + delta)*p_base, p_ceil) for delta in pref_profile.deltas])
                    if wf_state == self.WfState.BESS_DMAT_WTG_PZERO:
                        row_sect = {'Pref_Wind_MW_v': 0,
                                    'Pref_BESS_MW_v': pref_v,
                                    'Pref_BESS_MW_t': pref_profile.time_steps}
                    elif wf_state == self.WfState.WTG_DMAT_BESS_PMAX:
                        row_sect = {'Pref_Wind_MW_v': pref_v,
                                    'Pref_BESS_MW_v': self.p_max_bess,
                                    'Pref_Wind_MW_t': pref_profile.time_steps}
                    elif wf_state == self.WfState.WTG_DMAT_BESS_PZERO:
                        row_sect = {'Pref_Wind_MW_v': pref_v,
                                    'Pref_BESS_MW_v': 0,
                                    'Pref_Wind_MW_t': pref_profile.time_steps}
                    elif wf_state == self.WfState.WTG_DMAT_BESS_PMIN:
                        row_sect = {'Pref_Wind_MW_v': pref_v,
                                    'Pref_BESS_MW_v': -self.p_max_bess,
                                    'Pref_Wind_MW_t': pref_profile.time_steps}
                    elif wf_state == self.WfState.BESS_PNEG_DMAT_WTG_PZERO:
                        new_pref_v=[]
                        for i_pref_v in pref_v:
                            new_pref_v.append(-1*i_pref_v)
                        row_sect = {'Pref_Wind_MW_v': 0,
                                    'Pref_BESS_MW_v': new_pref_v,
                                    'Pref_BESS_MW_t': pref_profile.time_steps}
                    else:
                        raise CalcSheetError("wf_state")
                # add value
                else:
                    pref_v = min(row["Active Power (pu)"]*p_base, p_ceil)
                    if wf_state == self.WfState.BESS_DMAT_WTG_PZERO:
                        row_sect = {'Pref_Wind_MW_v': 0,
                                    'Pref_BESS_MW_v': pref_v}
                    elif wf_state == self.WfState.WTG_DMAT_BESS_PMAX:
                        row_sect = {'Pref_Wind_MW_v': pref_v,
                                    'Pref_BESS_MW_v': self.p_max_bess}
                    elif wf_state == self.WfState.WTG_DMAT_BESS_PZERO:
                        row_sect = {'Pref_Wind_MW_v': pref_v,
                                    'Pref_BESS_MW_v': 0}
                    elif wf_state == self.WfState.WTG_DMAT_BESS_PMIN:
                        row_sect = {'Pref_Wind_MW_v': pref_v,
                                    'Pref_BESS_MW_v': -self.p_max_bess}
                    elif wf_state == self.WfState.BESS_PNEG_DMAT_WTG_PZERO:
                        row_sect = {'Pref_Wind_MW_v': 0,
                                    'Pref_BESS_MW_v': -1 * pref_v}
                    else:
                        raise CalcSheetError("wf_state")
                row_sect.update({'DIR': f'{wf_state.name}_{new_test["DIR"]}'})
                row_sect.update({'File_Name': f'{wf_state.name}_{new_test["File_Name"]}'})
                row_sect_list.append(row_sect)
        else:
            raise CalcSheetError("p ref")
        return row_sect_list
    
    # if reactive power is specified, then add this in.
    def add_q_specs(self, row: pd.DataFrame, new_test: dict):
        row_sect = dict()
        # Reactive power is specified
        if 'Reactive Power (pu)' in row:
            row_sect = {'Init_Qpoc_pu_v': row["Reactive Power (pu)"]}
        else:
            row_sect = {'Init_Qpoc_pu_v': 0}
        return [row_sect]
            
    def add_v_specs(self, row: pd.DataFrame, new_test: dict):
        row_sect_list = []
        row_sect = dict()
        # Use data which has already been added to this test to calculate the values to be added
        # If the active power reference is a list, then use the first one as a reference
        if type(new_test["Pref_Wind_MW_v"]) == list:
            pref_wind = new_test["Pref_Wind_MW_v"][0]
        else:
            pref_wind = new_test["Pref_Wind_MW_v"]
        if type(new_test["Pref_BESS_MW_v"]) == list:
            pref_bess = new_test["Pref_BESS_MW_v"][0]
        else:
            pref_bess = new_test["Pref_BESS_MW_v"]
        ppoc = (pref_wind + pref_bess)/self.p_nom
        qpoc = new_test["Init_Qpoc_pu_v"]
        # if not an infinite bus, then we need to calculate the slack voltage from the POC voltage, otherwise it's the same
        if new_test["Infinite_Grid_v"] == 0:
            # if SCR1 FRT test, then use the pre fault SCR
            grid_fault_level = new_test["Grid_MVA_v"][0] if type(new_test["Grid_MVA_v"]) == list else new_test["Grid_MVA_v"]
            Zs = self.calc_grid_impedence(pu=True, fl=int(grid_fault_level), x2r=new_test["Grid_X2R_v"])
        elif new_test["Infinite_Grid_v"] == 1:
            Zs = None
        # General voltage profile specified
        if 'Vgrid Deltas (pu)' in row and 'Time Steps (s)' in row:
            init_vpoc = row["Voltage POC (pu)"] if 'Voltage POC (pu)' in row else self.v_set
            vgrid_profile = self.Profile(self.figure_references)
            vslack_starting_point = self.calc_vslack_from_vpoc(vpoc=init_vpoc, ppoc=ppoc, qpoc=qpoc, Zs=Zs)
            # split vgrid into a list of profiles
            vgrid_strs = row["Vgrid Deltas (pu)"].split("; ")
            timestep_strs = row["Time Steps (s)"].split("; ")
            if 'Vgrid Ramp (Hz/s)' in row: # add ramp information
                ramp_strs = row["Vgrid Ramp (Hz/s)"].split("; ")
                for index, vgrid_str in enumerate(vgrid_strs):
                    vgrid_profile.read_profile(v_data=vgrid_str,t_data=timestep_strs[index],r_data=ramp_strs[index])
                    if self.vgrid_change_method == "Vslack":
                        # vgrid change is applied via vslack
                        row_sect = {'Init_Vpoc_pu_v': self.calc_vpoc_from_vslack(vslack=vgrid_profile.deltas[0] + vslack_starting_point,ppoc=ppoc,qpoc=qpoc,Zs=Zs),
                                    'Vslack_pu_v': [delta + vslack_starting_point for delta in vgrid_profile.deltas],
                                    'Vslack_pu_t': vgrid_profile.time_steps,
                                    'Vslack_pu_r': vgrid_profile.ramps}
                    elif self.vgrid_change_method == "POC disturbance":
                        # vgrid change is applied at the POC via a dummy transformer
                        row_sect = {'Init_Vpoc_pu_v': init_vpoc,
                                    'Vslack_pu_v': vslack_starting_point,
                                    'Vpoc_disturbance_pu_v': vgrid_profile.deltas,
                                    'Vpoc_disturbance_pu_t': vgrid_profile.time_steps,
                                    'Vpoc_disturbance_pu_r': vgrid_profile.ramps}
                    else:
                        raise CalcSheetError("vgrid method")
                    row_sect_list.append(row_sect)
            else: # no ramp given
                for index, vgrid_str in enumerate(vgrid_strs):
                    vgrid_profile.read_profile(v_data=vgrid_str,t_data=timestep_strs[index])
                    if self.vgrid_change_method == "Vslack":
                        row_sect = {'Init_Vpoc_pu_v': self.calc_vpoc_from_vslack(vslack=vgrid_profile.deltas[0] + vslack_starting_point,ppoc=ppoc,qpoc=qpoc,Zs=Zs),
                                    'Vslack_pu_v': [delta + vslack_starting_point for delta in vgrid_profile.deltas],
                                    'Vslack_pu_t': vgrid_profile.time_steps}
                    elif self.vgrid_change_method == "POC disturbance":
                        row_sect = {'Init_Vpoc_pu_v': init_vpoc,
                                    'Vslack_pu_v': vslack_starting_point,
                                    'Vpoc_disturbance_pu_v': vgrid_profile.deltas,
                                    'Vpoc_disturbance_pu_t': vgrid_profile.time_steps}
                    else:
                        raise CalcSheetError("vgrid method")
                    row_sect_list.append(row_sect)
        elif 'Vgrid Deltas (pu)' in row: # if there are deltas but no time steps we need to fix the config sheet
            raise CalcSheetError("time steps")
        # There is a Vgrid profile in DMAT tests 149 - 166
        elif 'Event' in row and row["Event"] == "Vgrid":
            vgrid_profile = self.Profile(self.figure_references)
            vgrid_profile.read_profile(v_data=row["Delta (pu)"], t_data=row["Time Steps (s)"])
            # vslack starting point is calculated from Vpoc
            vslack_starting_point = self.calc_vslack_from_vpoc(vpoc=self.v_set, ppoc=ppoc, qpoc=qpoc, Zs=Zs)
            if self.vgrid_change_method == "Vslack":
                # vgrid change is applied via vslack
                row_sect = {'Init_Vpoc_pu_v': self.v_set,
                            'Vslack_pu_v': [delta + vslack_starting_point for delta in vgrid_profile.deltas],
                            'Vslack_pu_t': vgrid_profile.time_steps}
            elif self.vgrid_change_method == "POC disturbance":
                # vgrid change is applied at the POC via a dummy transformer
                row_sect = {'Init_Vpoc_pu_v': self.v_set,
                            'Vslack_pu_v': vslack_starting_point,
                            'Vpoc_disturbance_pu_v': vgrid_profile.deltas,
                            'Vpoc_disturbance_pu_t': vgrid_profile.time_steps}
            else:
                raise CalcSheetError("vgrid method")
            row_sect_list = [row_sect]
        # Vpoc is specified, calculate Vgrid/Vslack
        elif 'Voltage POC (pu)' in row:
            row_sect = {'Init_Vpoc_pu_v': row["Voltage POC (pu)"],
                        'Vslack_pu_v': self.calc_vslack_from_vpoc(vpoc=row["Voltage POC (pu)"], ppoc=ppoc, qpoc=qpoc, Zs=Zs)}
            row_sect_list = [row_sect]
        else:
            row_sect = {'Init_Vpoc_pu_v': self.v_set,
                        'Vslack_pu_v': self.calc_vslack_from_vpoc(vpoc=self.v_set, ppoc=ppoc, qpoc=qpoc, Zs=Zs)}
            row_sect_list = [row_sect]
        return row_sect_list
    
    def add_freq_specs(self, row: pd.DataFrame):
        row_sect = dict()
        row_sect_list = []
        if 'Freq Deltas (Hz)' in row:
            if 'Time Steps (s)' in row:
                freq_profile = self.Profile(self.figure_references)
                freq_strs = row["Freq Deltas (Hz)"].split("; ")
                time_strs = row["Time Steps (s)"].split("; ")
                if 'Freq Ramp (Hz/s)' in row:
                    ramp_strs = row["Freq Ramp (Hz/s)"].split("; ")
                    for index, freq_str in enumerate(freq_strs):
                        freq_profile.read_profile(v_data=freq_str,t_data=time_strs[index],r_data=ramp_strs[index])
                        freq_v  = [delta + self.f_nom for delta in freq_profile.deltas]
                        row_sect = {'Grid_Hz_t': freq_profile.time_steps,
                                    'Grid_Hz_v': freq_v,
                                    'Grid_Hz_r': freq_profile.ramps}
                        row_sect_list.append(row_sect)
                else:
                    for index, freq_str in enumerate(freq_strs):
                        freq_profile.read_profile(v_data=freq_str,t_data=time_strs[index])
                        freq_v  = [delta + self.f_nom for delta in freq_profile.deltas]
                        row_sect = {'Grid_Hz_t': freq_profile.time_steps,
                                    'Grid_Hz_v': freq_v}
                        row_sect_list.append(row_sect)
            else:
                raise CalcSheetError("time steps")
        elif 'Freq Ramp (Hz/s)' in row:
            raise CalcSheetError("freq deltas")
        else:
            row_sect_list = [{'Grid_Hz_v': self.f_nom}]
        return row_sect_list
    
    def add_qref_specs(self, row: pd.DataFrame, new_test: dict):
        row_sect = dict()
        if 'Event' in row and row["Event"] == "Qref":
            qref_profile = self.Profile(self.figure_references)
            qref_profile.read_profile(v_data=row["Delta (pu)"], t_data=row["Time Steps (s)"])
            q_starting_point = new_test["Init_Qpoc_pu_v"]
            row_sect = {'Qref_MVAr_v': [(delta + q_starting_point) * self.q_nom for delta in qref_profile.deltas],
                        'Qref_MVAr_t': qref_profile.time_steps}
        return [row_sect]
    
    def add_vref_specs(self, row: pd.DataFrame, new_test: list):
        row_sect = dict()
        if 'Event' in row and row["Event"] == "Vref":
            vref_profile = self.Profile(self.figure_references)
            vref_profile.read_profile(v_data=row["Delta (pu)"], t_data=row["Time Steps (s)"])
            # if the initial vref is not specified, then we use 1.034 pu
            v_starting_point = row["Vref (pu)"] if 'Vref (pu)' in row else self.v_set
            row_sect = {'Vref_pu_v': [delta + v_starting_point for delta in vref_profile.deltas],
                        'Vref_pu_t': vref_profile.time_steps}
            vref_init = v_starting_point + vref_profile.deltas[0]
        elif 'Vref Deltas (pu)' in row and 'Time Steps (s)':
            vref_profile = self.Profile(self.figure_references)
            # if the initial vref is not specified, then we use 1.034 pu
            v_starting_point = row["Vref (pu)"] if 'Vref (pu)' in row else self.v_set
            if 'Vref Ramp (pu/s)'  in row:
                vref_profile.read_profile(v_data=row["Vref Deltas (pu)"], t_data=row["Time Steps (s)"], r_data=row["Vref Ramp (pu/s)"])
                row_sect = {'Vref_pu_v': [delta + v_starting_point for delta in vref_profile.deltas],
                            'Vref_pu_t': vref_profile.time_steps,
                            'Vref_pu_r': vref_profile.ramps}
                vref_init = v_starting_point + vref_profile.deltas[0]
            else:
                vref_profile.read_profile(v_data=row["Vref Deltas (pu)"], t_data=row["Time Steps (s)"])
                row_sect = {'Vref_pu_v': [delta + v_starting_point for delta in vref_profile.deltas],
                            'Vref_pu_t': vref_profile.time_steps}
                vref_init = v_starting_point + vref_profile.deltas[0]
        # raise an error if there are vref deltas but no timesteps
        elif 'Vref Deltas (pu)' in row:
            raise CalcSheetError("time steps")
        # a single value may be given for vref
        elif 'Vref (pu)' in row:
            row_sect = {'Vref_pu_v': row["Vref (pu)"]}
            vref_init = row["Vref (pu)"]
        # if not otherwise specifed, then set the vref value depending on the values of Qpoc and Vpoc
        else:
            row_sect = {'Vref_pu_v': self.calc_vref_from_qpoc_and_vpoc(qpoc=new_test["Init_Qpoc_pu_v"],vpoc=new_test["Init_Vpoc_pu_v"])}
        
        #whether or not any deltas are applied, the following are only carried out if Vref is specified
        if 'Vref (pu)' in row:
            if 'Voltage POC (pu)' in row and 'Reactive Power (pu)' in row:
                raise CalcSheetError("vref, vpoc, and qpoc")
            # if Qpoc and Vref are specified but not Vpoc, then we reset Vpoc to Vpoc(Qpoc, Vref)
            elif 'Reactive Power (pu)' in row:
                row_sect.update({'Init_Vpoc_pu_v': self.calc_vpoc_from_qpoc_and_vref(qpoc=new_test["Init_Qpoc_pu_v"],vref=vref_init)})
            # if Vref and Vpoc are specified but Qpoc isn't then we reset Qpoc to Qpoc(Vref, Vpoc)
            elif 'Voltage POC (pu)' in row:
                row_sect.update({'Init_Qpoc_pu_v': self.calc_qpoc_from_vref_and_vpoc(vref=vref_init,vpoc=new_test["Init_Vpoc_pu_v"])})
        return [row_sect]
    
    def add_osc_specs(self, row: pd.DataFrame, new_test: dict):
        row_sect_list = []
        row_sect = dict()
        if 'Time Steps (s)' in row and 'Osc Freqs (Hz)' in row and 'Osc Magnitude (pu)' in row and 'Osc Phase (deg)' in row:
            # read the time steps
            try:
                time_steps = json.loads(row["Time Steps (s)"])
            except:
                raise CalcSheetError("time steps")
            timing_sig = [0, time_steps[0], time_steps[1], row["End Run (s)"]]
            # split the osc freqs into different test sets
            osc_freq_strs = str(row["Osc Freqs (Hz)"]).split("; ")
            # split the freq steps into different test sets
            steps = str(row["Step (Hz)"]).split("; ") if 'Step (Hz)' in row else None
            # iterrate through each test group
            for index, osc_freq_str in enumerate(osc_freq_strs):
                if not (steps is None or steps[index] == None):
                    osc_freq_steps = osc_freq_str.split(":")
                    if not len(osc_freq_steps) == 2:
                        raise CalcSheetError("osc freq steps")
                    step = float(steps[index])
                    # make a list of all the frequencies which should be added
                    osc_freq_range = []
                    temp = float(osc_freq_steps[0])
                    while temp < float(osc_freq_steps[1]):
                        osc_freq_range.append(temp)
                        temp = temp + step
                else: # no step specified
                    osc_freq_range = [osc_freq_str]
                # iterate through each frequency level and add a test
                for osc_freq_step in osc_freq_range:
                    row_sect = {"Vslack_osc_amplitude_v": [0, row["Osc Magnitude (pu)"], 0, 0],
                                "Vslack_osc_amplitude_t": timing_sig,
                                "Vslack_osc_Hz_v": [0, osc_freq_step, 0, 0],
                                "Vslack_osc_Hz_t": timing_sig,
                                "Vslack_osc_phase_deg_v": [0, row["Osc Phase (deg)"], 0, 0],
                                "Vslack_osc_phase_deg_t": timing_sig}
                    row_sect_list.append(row_sect)
        elif 'Osc Magnitude (pu)' in row or 'Osc Phase (deg)' in row or 'Osc Freqs (Hz)' in row:
            raise CalcSheetError("osc info")
        return row_sect_list
    
    def add_phase_specs(self, row: pd.DataFrame, new_test: dict):
        row_sect = dict()
        row_sect_list = []
        if 'Angle Change (deg)' in row and 'Apply Event (s)' in row:
            # Create a list of each phase which must be applied in a different test
            phase_strs = str(row["Angle Change (deg)"]).split("; ")
            for phase_str in phase_strs:
                row_sect = {"Grid_phase_degs_v": [0, float(phase_str), float(phase_str)],
                            "Grid_phase_degs_t": [0, row["Apply Event (s)"], row["End Run (s)"]]}
                row_sect_list.append(row_sect)
        elif 'Angle Change (deg)' in row:
            raise CalcSheetError("apply event")
        return row_sect_list
            
    
    def is_number(self, string:str):
        try:
            float(string)
            return True
        except ValueError:
            return False
    
    def add_fault_specs(self, row: pd.DataFrame, new_test: dict):
        row_sect = dict()
        row_sect_list = []
        if 'Apply Fault (s)' in row and 'Fault Duration (s)' in row:
            fault_duration = str(row["Fault Duration (s)"])
            if self.is_number(fault_duration):
                # In the SCR1 FRT case, the fault level changes, the fault impedance should be based on the fault level before the fault
                grid_fault_level = new_test["Grid_MVA_v"][0] if type(new_test["Grid_MVA_v"])== list else new_test["Grid_MVA_v"]
                grid_impedance = self.calc_grid_impedence(pu=1, fl=grid_fault_level, x2r=new_test["Grid_X2R_v"])
                # if the fault x2r is specified then we use this, otherwise we assume the grid x2r
                fault_x2r = row["Fault X/R"] if 'Fault X/R' in row else new_test["Grid_X2R_v"]
                # If both fault impedance and fault voltage are specified then we use fault voltage
                if 'Applied Fault Voltage (pu)' in row:
                    # if fault type is not specified then we assume 3 phase to ground
                    fault_type = row["Fault Type"] if 'Fault Type' in row else "3PHG"
                    fault_voltage_strs = str(row["Applied Fault Voltage (pu)"]).split("; ")
                    steps = str(row["Voltage Step (pu)"]).split("; ") if 'Voltage Step (pu)' in row else None
                    for index, fault_voltage_str in enumerate(fault_voltage_strs):
                        if not (steps is None or steps[index] == None): # if there are steps then we create a range and iterate through
                            fault_voltage_steps = fault_voltage_str.split(":")
                            if not len(fault_voltage_steps) == 2:
                                raise CalcSheetError("fault voltage steps")
                            step = float(steps[index])
                            # make a list of all the voltage steps which should be added
                            fault_voltage_range = []
                            temp = float(fault_voltage_steps[0])
                            while temp < float(fault_voltage_steps[1]):
                                fault_voltage_range.append(temp)
                                temp = temp + step
                        else: # when there are no steps the range is just the fault voltage given
                            fault_voltage_range = [fault_voltage_str]
                        # iterate through each step and add info
                        for fault_voltage_step in fault_voltage_range:
                            fault_impedance = self.calc_fault_impedance(fault_voltage=fault_voltage_step, fault_distance=1, grid_impedance=grid_impedance)
                            fault_impedance = fault_impedance/self.z_base
                            row_sect = {'Enable_Fault_Studies_v': 1,
                                        'Fault_Strategy_v': 1,
                                        'Fault_X2R_v': fault_x2r,
                                        'Fault_Type_v': self.get_fault_code(fault_type),
                                        'Fault_Timing_Signal_t': f'[0, {row["Apply Fault (s)"]}, {row["Apply Fault (s)"] + row["Fault Duration (s)"]}]',
                                        'Fault_Timing_Signal_v': f'[0, 1, 0]',
                                        'Zf2Zs_v': self.calc_fault_multiplier_from_fault_impedance(fault_impedance=fault_impedance,grid_impedance=grid_impedance)}
                            row_sect_list.append(row_sect)
                elif 'Fault Impedance (pu)' in row:
                    # Temporary over voltage specfied as "Yf = jXc (U_Ov = 1.15pu)"
                    if "Yf = jXc" in row["Fault Impedance (pu)"] and "U_Ov =" in row["Fault Impedance (pu)"]:
                        # 'Yf = jXc (U_Ov = 1.15pu)' --> ['Yf', 'jXc (U_Ov', '1.15pu') ]
                        fault_voltage = row["Fault Impedance (pu)"].split("=")
                        # '1.15pu)' --> ['1.15', ')']
                        fault_voltage = fault_voltage[2].split("pu")[0]
                        fault_impedance = self.calc_fault_impedance(fault_voltage=float(fault_voltage), fault_distance=1, grid_impedance=grid_impedance)
                        row_sect = {'TOV_Shunt_uF_v': 1/(fault_impedance * self.f_nom * 2 * math.pi) * 1e6,
                                    'TOV_Timing_Signal_t': f'[0, {row["Apply Fault (s)"]}, {row["Apply Fault (s)"] + row["Fault Duration (s)"]}]',
                                    'TOV_Timing_Signal_v': f'[0, 1, 0]'}
                        row_sect_list = [row_sect]
                    else:
                        # if fault type is not specified then we assume 3 phase to ground
                        fault_type = row["Fault Type"] if 'Fault Type' in row else "3PHG"
                        # Enable fault studies should be set to one to allow fault to be added
                        # Fault strategy=1 means that the "Zs2Zf method" is used which means that the impedance is calculated as a ratio of the grid impedance.
                        # The alternative is fault startegy =0 which means the fault impedance is calculated based off Udip in PSCAD.
                        row_sect = {'Enable_Fault_Studies_v': 1,
                                    'Fault_Strategy_v': 1,
                                    'Fault_X2R_v': fault_x2r,
                                    'Fault_Type_v': self.get_fault_code(fault_type),
                                    'Fault_Timing_Signal_t': f'[0, {row["Apply Fault (s)"]}, {row["Apply Fault (s)"] + row["Fault Duration (s)"]}]',
                                    'Fault_Timing_Signal_v': f'[0, 1, 0]'}
                        # Resistive faults
                        if row["Fault Impedance (pu)"] == "Zf=Rf = 1, 5 and 10 Ohm":
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
                            row_sect.update({'Zf2Zs_v': self.read_fault_multiplier(row["Fault Impedance (pu)"], Zs=grid_impedance)})
                            row_sect_list.append(row_sect)
                elif self.is_list(new_test["Grid_MVA_v"]):
                    pass
                else: # neither fault impedance or fault voltage are specified.
                    ic(new_test["Grid_MVA_v"])
                    raise CalcSheetError("fault impedance or fault voltage")
            # multiple fault ride through
            elif fault_duration in ["S1", "S2", "S3", "S4", "S5"]:
                row_sect = {'Enable_Fault_Studies_v': 1,
                            'Fault_Strategy_v': 0,}
                row_sect.update(self.mfrt_seq(row, new_test))
                row_sect_list.append(row_sect)
            else: # fault duration is not a number and not  and mfrt code
                raise CalcSheetError("fault duration")
        elif 'Apply Fault (s)' in row:
            raise CalcSheetError("fault duration")
        elif 'Fault Duration (s)' in row:
            raise CalcSheetError("apply fault")
        return row_sect_list
        
##################### MFRT GENERATOR #######################

    def mfrt_seq(self, row: pd.DataFrame, new_test: dict):
        random.seed(row["Test Number"])
        fault_info = []
        fault_type_options = ["3PHG", "2PHG", "1PHG", "L-L"]
        if row["Fault Duration (s)"] == "S1":
            fault_info = [("3PHG", 5, 0.1, 0.25),
                        ("3PHG", 5.25, 0.1, 0.25),
                        ("3PHG", 5.5, 0.1, 0.25),
                        ("2PHG", 8, 0.1, 0.25),
                        ("2PHG", 11, 0.1, 0.25),
                        ("2PHG", 13, 0.1, 0.25)]
            time = [13]
            fault_duration = 0.1
        else:
            # From S5.2.5.5:
            # Up to 15 faults in any 5 minute period (arbitrarily I decided to constrain it to between 5 and 15 faults in 2 to 5 mins.)
            # up to 6 disturbances cause the POC voltage to drop below 50%
            # minimum clearance from the end of one disturbance and commencement of the next disturbance may be 0 ms
            # cumulative time that voltage at the connection point is lower than 90% of normal voltage will not exceed 1800 ms
            # OR the time integral of difference between V dip and normal V when V dip > %90* Vnom will not exceed 1s
            min_no_faults = 5
            max_no_faults = 15
            max_time = 300
            # Choose a ranom number of faults to apply between 5 and 15
            no_faults = random.randint(min_no_faults, max_no_faults)
            # We can specify a minimum time between faults
            min_time_between_faults = 0.001
            # By using the integer random sampling method we can control the level of precision considered.
            accuracy = 0.001
            time_range = range(0, int(min(max_time, row["End Run (s)"])/accuracy), int(min_time_between_faults/accuracy))
            time = random.sample(time_range, k=no_faults)
            time.sort()
            # Convert the list back to seconds.
            time = [item*accuracy for item in time]
            # Calculate the fault multiplier for which the POC voltage is 50%
            # TODO
            no_vdips_gt50 = 0
            # Calculate the fault multiplier for which the POC voltage is 90%
            # TODO
            time_vdips_gt_10 = 0
            # Loop through each time step.
            for i in range(len(time)):
                fault_type = random.choice(fault_type_options)
                next_time = int(row["End Run (s)"]) if i == len(time)-1 else time[i+1]
                fault_duration = random.randrange(1, int((next_time-time[i])/accuracy), 1)*accuracy
                fault_multiplier = random.randint(0, 5)
                fault_info.append((fault_type, time[i], fault_duration, fault_multiplier))
        row_sect = self.generate_mfrt(fault_info)
        row_sect.update({'Fault_X2R_v': new_test["Grid_X2R_v"]})
        return row_sect
    
    def generate_mfrt(self, fault_info: list):
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
            if index > 0 and index < len(fault_info) - 1:
                post_fault_durations.append(fault_time - fault_times[index-1] - fault_durations[index-1])
            fault_zs_multiplier.append(fault_impedence_multiplier)
            fault_timing_signal_t.extend([fault_time, fault_time+fault_duration])
            fault_timing_signal_v.extend([1, 0])
            fault_times.append(fault_time)
        
        row_sect = {'Fault_Timing_Signal_t': fault_timing_signal_t,
                    'Fault_Timing_Signal_v': fault_timing_signal_v,
                    'Zf2Zs_t': fault_times,
                    'Zf2Zs_v': fault_zs_multiplier,
                    'Fault_Type_t': fault_times,
                    'Fault_Type_v': fault_types}    
        return row_sect
    
    def calc_ramp(self, row: pd.DataFrame):
        d_freq = abs(row["Freq Target"] - self.f_nom)
        d_t = d_freq/row["Freq Tamp Hz/s"]
        t_end_ramp = min(row["Apply Step (s)"] + d_t, row["End Run (s)"])
        return t_end_ramp 
        
##################### FAULT CALCULATIONS #######################
    def calc_fault_multiplier_from_fault_impedance(self, fault_impedance: float, grid_impedance: tuple):
        (_, _, Zs) = grid_impedance
        fault_multiplier = fault_impedance/Zs
        return fault_multiplier
        
    def calc_fault_impedance(self, fault_voltage: float, fault_distance: float, grid_impedance: tuple):
        (_, _, grid_impedance) = grid_impedance
        grid_impedance = grid_impedance*self.z_base
        fault_impedance = fault_distance*grid_impedance*fault_voltage/abs(1 - fault_voltage)
        return fault_impedance
    
    def read_fault_multiplier(self, zf_str: str, Zs: tuple):
        # can't set fault impedance to 0 so set it to a very low number.
        min_fault_impedance = 0.000000001
        str_list = zf_str.split("=")
        # If there is a "[comment]" part, remove this.
        if "[" in str_list[1]:
            str_list[1] = str_list[1].split("[")[0]
            str_list = str_list[0:2]
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
                if len(str_list) == 2 and str_list[1] == "Zs" and self.is_number(str_list[0]):
                    # "2xZs" str_list = [2, Zs]
                    multiplier = float(str_list[0])
                elif self.is_number(str_list[1]):
                    Zf = float(str_list[1])
                    multiplier = self.calc_fault_multiplier_from_fault_impedance(fault_impedance=Zf, grid_impedance=Zs)
                else:
                    raise CalcSheetError("fault impedance")
        else:
            raise CalcSheetError("fault impedance")
        return multiplier
    
    def calc_fault_voltage(self, Zf: tuple, Zs: tuple):
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
    def is_list(self, string: str):
        try:
            _ = json.loads(str(string))
            return True
        except:
            return False
            
    def read_scr_and_x2r(self, row_scr: str, row_x2r: str):
        withstand_scr = self.system_inf.loc["Withstand SCR"]["Var Val"]
        # check withstand scr is a number
        if not self.is_number(withstand_scr):
            raise CalcSheetError("withstand scr")
        row_scr = str(row_scr).replace("WITHSTAND SCR", str(withstand_scr))
        scr_strs = str(row_scr).split("; ")
        x2r_strs = str(row_x2r).split("; ")
        scr_and_x2r = []
        poc_scr = json.loads(self.system_inf.loc["POC SCR"]["Var Val"])
        poc_x2r = json.loads(self.system_inf.loc["POC XR ratio"]["Var Val"])
        
        for scr_str in scr_strs:
            for x2r_str in x2r_strs:
                # get a list of x2r values
                x2rs = []
                if x2r_str.isnumeric():
                    x2rs = [float(x2r_str)]
                elif x2r_str == "POC":
                    x2rs = [float(poc_x2r[0]), float(poc_x2r[1])]
                elif x2r_str == "POC MIN":
                    x2rs = [float(poc_x2r[0])]
                elif x2r_str == "POC MAX":
                    x2rs = [float(poc_x2r[1])]
                else:
                    raise CalcSheetError("X2R")
                # add x2r info to scr info
                if self.is_number(scr_str):
                    for x2r in x2rs:
                        scr_and_x2r.append((float(scr_str), x2r))
                elif scr_str == "POC":
                    # if both SCR and X/R are POC, then we only add two rows
                    if x2r_str == "POC":
                        scr_and_x2r.append((float(poc_scr[0]), x2rs[0]))
                        scr_and_x2r.append((float(poc_scr[1]), x2rs[1]))
                    else:
                        for x2r in x2rs:
                            scr_and_x2r.append((float(poc_scr[0]), x2r))
                            scr_and_x2r.append((float(poc_scr[1]), x2r))
                elif scr_str == "POC MIN":
                    for x2r in x2rs:
                        scr_and_x2r.append((float(poc_scr[0]), x2r))
                elif scr_str == "POC MAX":
                    for x2r in x2rs:
                        scr_and_x2r.append((float(poc_scr[1]), x2r))
                elif self.is_list(scr_str):
                    for x2r in x2rs:
                        scr_and_x2r.append((scr_str, x2r))
                else:
                    raise CalcSheetError("SCR")
        return scr_and_x2r
    
    def calc_fault_level(self, scr: float):
        fault_level = scr * self.p_nom
        return fault_level
    
    def calc_grid_impedence(self, pu: bool, fl, x2r: float):
        z = self.v_nom**2 / fl
        r = z / math.sqrt(1+x2r**2)
        x = r * x2r
        if pu == True:
            r = r/self.z_base
            x = x/self.z_base
            z = z/self.z_base
        l = x / (2 * self.f_nom * math.pi)
        return (r,x,z)
        
    def calc_vslack_from_vpoc(self, vpoc: float, ppoc: float, qpoc: float, Zs: tuple):
        # if infinite grid, then vlsack = vpoc
        if Zs == None:
            return vpoc
        spoc = math.sqrt(ppoc**2 + qpoc**2)
        (rgrid_pu, xgrid_pu, zgrid_pu) = Zs

        vslack = math.sqrt(vpoc**2 + spoc**2/vpoc**2*zgrid_pu**2 - 2*(ppoc*rgrid_pu + qpoc*xgrid_pu))
        return vslack
        
    def calc_vpoc_from_vslack(self, vslack: float, ppoc: float, qpoc: float, Zs: tuple):
        # if infinite grid, then vlsack = vpoc
        if Zs == None:
            return vslack
        spoc = math.sqrt(ppoc**2 + qpoc**2)
        (rgrid_pu, xgrid_pu, zgrid_pu) = Zs
        vgrid = math.sqrt((vslack**2 + 2*(ppoc*rgrid_pu + qpoc*xgrid_pu) + math.sqrt((vslack**2 + 2*(ppoc*rgrid_pu + qpoc*xgrid_pu))**2 - 4*spoc**2*zgrid_pu**2))/2)
        return vgrid
    
    def calc_qpoc_from_vref_and_vpoc(self, vref: float, vpoc: float):
        qpoc = (vref - vpoc)/self.qv_droop
        return qpoc
    
    def calc_vref_from_qpoc_and_vpoc(self, qpoc: float, vpoc: float):
        vref = vpoc + qpoc*self.qv_droop
        return vref
    
    def calc_vpoc_from_qpoc_and_vref(self, qpoc: float, vref: float):
        vpoc = -qpoc*self.qv_droop + vref
        return vpoc


    # To simplify access to a profile, there is a profile subclass which holds the information for a signal including _v, _t, and _r information
    class Profile(object):
        def __init__(self, figure_references: pd.DataFrame):
            self.figure_references = figure_references
            self.deltas = None
            self.time_steps = None
            self.ramps = None
            self.profile_type = None

        def check_valid_input(self, v_data:str, t_data:str, r_data=None):
            # check if the data inputs can be read as a list
            try:
                _ = json.loads(v_data)
                return "manual"
            except json.decoder.JSONDecodeError:
                ic(v_data)
                # check that the data inputs are all the same
                if not v_data == t_data:
                    ic()
                    return "invalid"
                if not r_data is None:
                    ic(r_data)
                    if not v_data == r_data:
                        return "invalid"
                # check if the data inputs can be read as a default figure
                try:
                    ic()
                    self.get_default_profile_type(v_data)
                    return "default"
                except KeyError:
                    return "invalid"
        
        def get_default_profile_type(self, figure_name: str):
            ic(figure_name)
            figure_reference = self.figure_references.loc[figure_name]
            self.profile_type =  figure_reference["Type"]
            
        def read_default_profile(self, fig: str, ramp = False):
            if not self.check_valid_input(v_data=fig, t_data=fig) == "default":
                raise CalcSheetError("figure name")
            # default profile is refered to as "fig"
            figure_reference = self.figure_references.loc[fig]
            if not self.check_valid_input(v_data=figure_reference["Deltas"], t_data=figure_reference["Time Steps"]) == "manual":
                raise CalcSheetError("figure reference")
            self.deltas = json.loads(figure_reference["Deltas"])
            self.time_steps = json.loads(figure_reference["Time Steps"])
            if ramp:
                self.ramps = figure_reference["Ramp"] # we can just read the ramp in as a string
        
        def read_profile(self, v_data: str, t_data: str, r_data=None):
            v_data = str(v_data).replace(" ","")
            t_data = str(t_data).replace(" ","")
            input_type = self.check_valid_input(v_data=v_data, t_data=t_data, r_data=r_data)
            # manual profile is refered to as [1,2,3]
            if input_type == "manual":
                self.deltas = json.loads(v_data)
                self.time_steps = json.loads(t_data)
                if not r_data is None:
                    self.ramps = str(r_data).replace(" ","") # we can just log ramp data as a string
            elif input_type == "default":
            # if the profile can't be recognised, then try read the default profile
            # If there is a ramp in the profile, then read the ramp of the default profile
                if r_data is None:
                    self.read_default_profile(fig=v_data, ramp=True)
                else:
                    self.read_default_profile(fig=v_data, ramp=True)
            elif input_type == "invalid":
                raise CalcSheetError("time step, ramp, or delta")
            
