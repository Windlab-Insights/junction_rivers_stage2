import pandas as pd
import json
from typing import Dict
import math
from icecream import ic


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
                        new_row = self.add_new_row(row)
                        ic(new_row)
                        new_row.update({'File_Name': f'{category_formated}_test_{row["Test"].replace("-","_")}'})
                        new_row.update({'DIR': category_formated})
                        new_row.update({'En_SMIB_init_v': 1})
                        new_row.update({'Infinite_Grid_v': 1})
                        ic(new_row)
                        spec_df = spec_df.append(new_row, ignore_index=True)
                        ic(spec_df)
                        
        spec_df.to_csv(spec_path)
        
        
################# ADD PARAMS FOR A TEST ################

    def add_new_row(self, row):
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
        
        return new_row
    
    def add_file_info(self, row):
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
    
    def add_v_specs(self, row):
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
            # TODO
            print(f"## Throw an exception")
        ic(row_sect)
        return row_sect
    
    def add_q_specs(self, row):
        ic('add_q_specs')
        row_sect = dict()
        row_sect = {'Init_Qpoc_pu_v': self.q_set}
        ic(row_sect)
        return row_sect
    
    def add_p_specs(self, row, wf_state):
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
    
    def add_freq_specs(self, row):
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
    
    def add_fault_specs(self, row):
        ic("add_fault_specs")
        row_sect = dict()
        if 'Fault type' in row:
            if 'Fault Impedance' in row:
                if 'Fault duration(s)' in row:
                    if 'Apply Fault(s)' in row:
                        row_sect = {'Init_Fault_MVA': self.calc_fault_level(row),
                                    'Init_Fault_X_on_R': self.read_x2r(row),
                                    'Fault_Type_v': self.get_fault_code(row["Fault type"]),
                                    'Fault_Depth': 0,
                                    'Fault_Time': row["Apply Fault(s)"],
                                    'Fault_Duration': row["Fault duration(s)"],
                                    'Post_Fault_Duration': row["End Run (s)"] - row["Apply Fault(s)"] - row["Fault duration(s)"],
                                    'Enable_Fault_Studies_v': 1,
                                    'Fault_Strategy_v': 0,
                                    'Fault_X2R_v': self.read_x2r(row),
                                    'Fault_Timing_Signal_t': f'[0, {row["Apply Fault(s)"]}, {row["Apply Fault(s)"] + row["Fault duration(s)"]}',
                                    'Fault_Timing_Signal_v': f'[0, 1, 0]',
                                    'Ures_v': 0}
                        
        return row_sect
        
##################### MORE TOOLS #######################

    def generate_seq1(self,row):
        row_sect = dict()
        row_sect = {'Fault_Types': [self.get_fault_code(["3PHG", "3PHG", "3PHG", "2PHG", "2PHG", "2PHG"])],
                    'Fault_Durations': [0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
                    'Post_Fault_Durations': [0.15],
                    'Fault_Zs_Multiplier': ,
                    'Fault_Timing_Signal_t': ,
                    'Fault_Timing_Signal_v': [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
                    'Zf2Zs_v': ,
                    'Zf2Zs_t': ,
                    'Fault_ Type_t': ,
                    'Fault_Type_v': }
        return row_sect
    
    def calc_ramp(self, row):
        ic('calc_ramp')
        d_freq = abs(row["Freq target"] - self.f_nom)
        d_t = d_freq/row["Freq ramp Hz/s"]
        t_end_ramp = min(row["Apply step (s)"] + d_t, row["End Run (s)"])
        return t_end_ramp
        
        
##################### FAULT CALCULATIONS #######################
    def calc_fault_impedance(self, row):
        if row["Fault Impedance"] == "Zf=0":
            fault_impedance = 0
        elif row["Fault Impedance"] == "Zf=Zs":
            fault_impedance = self.calc_fault_impedance(row, 1)
        return fault_impedance
    
    def calc_fault_voltage(self,row):
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
    def get_fault_code(self, code):
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
    def read_scr(self, row):
        ic("read_scr")
        ic(row)
        scr = row["SCR"]
        ic(scr)
        if scr.isnumeric():
            return scr
        elif scr == "POC":
            ic(json.loads(self.system_inf.loc["POC SCR"]["Var Val"])[1])
            return json.loads(self.system_inf.loc["POC SCR"]["Var Val"])[1]
        else:
            # TODO add exception
            ic("error")
            return json.loads(self.system_inf.loc["POC SCR"]["Var Val"])[1]
        
    def read_x2r(self, row):
        ic("read_scr")
        x2r = row["X/R"]
        ic(x2r)
        if x2r.isnumeric():
            return x2r
        elif x2r == "POC":
            ic(self.system_inf.loc["POC XR ratio"]["Var Val"])
            return self.system_inf.loc["POC XR ratio"]["Var Val"]
        else:
            # TODO add exception
            ic("error")
            return self.system_inf.loc["POC XR ratio"]["Var Val"]
    
    def calc_fault_level(self, row):
        ic("calc_fault_level")
        scr = self.read_scr(row)
        fault_level = scr * self.p_nom
        ic(fault_level)
        return fault_level
    
    def calc_grid_impedence(self, row, pu):
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
        
    def calc_vslack_from_vpoc(self, row):
        ic("calc_vslack_from_vpoc")
        vpoc = row["Vpoc"]
        ppoc = row["Active Power (pu)"]
        qpoc = self.q_set
        spoc = math.sqrt(ppoc**2 + qpoc**2)
        (rgrid_pu, xgrid_pu, zgrid_pu) = self.calc_grid_impedence(row, 1)
        vslack = math.sqrt(vpoc**2 + spoc**2/vpoc**2*zgrid_pu**2 - 2*(ppoc*rgrid_pu + qpoc*xgrid_pu))
        return vslack
        
    def calc_vpoc_from_vslack(self, row):
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
    
    