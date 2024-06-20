import pandas as pd
import json
from typing import Dict


class SpecGenerator():
    
    HEADERS = [
        'Grouping',
        'Key_Tests',
        'File_Name',
        'Post_Init_Duration_s',
        'Category',
        'DIR',
        ]
    
    
    def __init__(self):
        self.start = 1

        
    def spec_generator (self, calc_sheet_path, spec_path, mapping_sheet):
        
        checklist = pd.read_excel(calc_sheet_path, sheet_name="Checklist", index_col="Checklist", header=0)
        print(f"### checklist = {checklist}")
        print(f"### headers = {checklist.columns.values}")
        self.system_inf = pd.read_excel(calc_sheet_path, sheet_name="System_inf", index_col="Var name", header=0)
        with open(mapping_sheet, 'r') as f:
            sim_headers = list(json.load(f).keys())
            f.close()
        print(sim_headers)
        spec_df = pd.DataFrame()
        print(f"### checklist flat run = {checklist.loc['Flat Run']}")
        print(f"### flat run = {checklist.loc['Flat Run']['Action']}")
        
        if checklist.loc["Flat Run"]["Action"] == "Yes":
            
            print(f"### in if statement")
            flat_run = pd.read_excel(calc_sheet_path, sheet_name="Flat Run", header=0)
            print(f"### flat run = {flat_run}")
            print(f"### flat run headers = {flat_run.columns.values}")
            
            for index, row in flat_run.iterrows():
                print(f"### row = {row}")
                print(f"### row action {row['Action']}")
                print(f"### spec df 0 : {spec_df}")  
                if row['Action'] == "Yes":
                    new_row ={
                        'Grouping': '',
                        'Key_Tests': '',
                        'File_Name': f'flat_run_test',
                        'Post_Init_Duration_s': row['End Run'],
                        'Category': '',
                        'DIR': 'flat_run',
                        'Grid_SCR_v': row['SCR'],
                        'Grid_MVA_v': self.calc_fault_level(row['SCR']),
                        'Grid_X2R_v': row['X/R'],
                        'Init_Vpoc_pu_v': 1.034,
                        
                    }
                    spec_df = spec_df.append(new_row, ignore_index=True)
                    print(f"### spec df : {spec_df}")  
                    
        print(f"### spec df to : {spec_df}")           
        spec_df.to_csv(spec_path)
        print(f"### csv saved to : {spec_path}")
        
    def calc_fault_level(self, scr):
        
        p_max_poc = self.system_inf.loc['Nominal Power']['Var Val']
        fault_level = scr * p_max_poc
        
        return fault_level