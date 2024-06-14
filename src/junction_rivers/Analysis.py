import pandas as pd
from collections import OrderedDict
from typing import Union
from pathlib import Path
import matplotlib.pyplot as plt

NUM_WTG1= 11
NUM_WTG2= 14

class Analysis:
    
    """
    Freqency Disturbance Analysis
    """

    def fdroop (self, tuning, scenario: pd.DataFrame, run_data: pd.DataFrame):
        init_ppoc_wind_mw = float(scenario['Pref_Wind_MW_v'])
        init_ppoc_bess_mw = float(scenario['Pref_BESS_MW_v'])

        model_init_time = float(tuning['TIME_Full_Init_Time_sec'])
        post_init_time = model_init_time + float(scenario['Post_Init_Duration_s.1'])

        data = run_data[model_init_time:post_init_time]

        ppoc_target_mw = init_ppoc_wind_mw + init_ppoc_bess_mw

        analysis_map = OrderedDict()

        freq_during_disturbance_hz = data["Freq_POC"].values[-1]
        poc_p_during_disturbance_mw = data["P_POC"].values[-1]
        wtg_p_during_disturbance_mw = data["P_POI_WT"].values[-1]
        bess_p_during_disturbance_mw = data["P_POI_BESS"].values[-1]

        analysis_map['Freq_During_Dist_Hz'] = freq_during_disturbance_hz
        analysis_map['POC_P_During_Dist_MW'] = poc_p_during_disturbance_mw
        analysis_map['WTG_P_During_Dist_MW'] = wtg_p_during_disturbance_mw
        analysis_map['BESS_P_During_Dist_MW'] = bess_p_during_disturbance_mw

        analysis_map['Delta_Freq_During_Dist_Hz'] = (freq_during_disturbance_hz - 50)/50
        analysis_map['Delta_POC_P_During_Dist_MW'] = (poc_p_during_disturbance_mw - ppoc_target_mw)/150
        analysis_map['Delta_WTG_P_During_Dist_MW'] = (wtg_p_during_disturbance_mw - init_ppoc_wind_mw)/150
        analysis_map['Delta_BESS_P_During_Dist_MW'] = (bess_p_during_disturbance_mw - init_ppoc_wind_mw)/150
        
        self.fdroop_log = self.fdroop_log.append(other=pd.DataFrame([analysis_map]), ignore_index=True)

        return

    def plot_fdroop(self, pdf_path: Union[Path, str]):
        plt.clf()
        fig = plt.figure()
        ax_total = plt.subplot(2,1,1)
        ax_delta = plt.subplot(2,1,2)
        ax_total.set_title('Droop response')
        ax_total.set_xlabel('Frequency (Hz)')
        ax_total.set_ylabel('Active Power (MW)')
        ax_delta.set_title('Droop response in per unit')
        ax_delta.set_xlabel('Delta Frequency (pu)')
        ax_delta.set_ylabel('Delta Active Power (pu)')
        
        freq_during_disturbance_hz = self.fdroop_log['Freq_During_Dist_Hz']
        poc_p_during_disturbance_mw = self.fdroop_log['POC_P_During_Dist_MW']
        wtg_p_during_disturbance_mw = self.fdroop_log['WTG_P_During_Dist_MW']
        bess_p_during_disturbance_mw = self.fdroop_log['BESS_P_During_Dist_MW']
        
        delta_freq_during_disturbance_hz = self.fdroop_log['Delta_Freq_During_Dist_Hz']
        delta_poc_p_during_disturbance_mw = self.fdroop_log['Delta_POC_P_During_Dist_MW']
        delta_wtg_p_during_disturbance_mw = self.fdroop_log['Delta_WTG_P_During_Dist_MW']
        delta_bess_p_during_disturbance_mw = self.fdroop_log['Delta_BESS_P_During_Dist_MW']
        
        ax_total.plot(freq_during_disturbance_hz, poc_p_during_disturbance_mw, 'r', ls='', marker='o', label='POC')
        ax_total.plot(freq_during_disturbance_hz, wtg_p_during_disturbance_mw, 'b', ls='', marker='o', label='WT')
        ax_total.plot(freq_during_disturbance_hz, bess_p_during_disturbance_mw, 'g', ls='', marker='o', label='BESS')
        ax_total.legend(loc="upper left")
        
        ax_delta.plot(delta_freq_during_disturbance_hz, delta_poc_p_during_disturbance_mw,'r', ls='', marker='o', label='POC')
        ax_delta.plot(delta_freq_during_disturbance_hz, delta_wtg_p_during_disturbance_mw, 'b', ls='', marker='o', label='WT')
        ax_delta.plot(delta_freq_during_disturbance_hz, delta_bess_p_during_disturbance_mw, 'g', ls='', marker='o', label='BESS')
        ax_total.legend(loc="upper left")
        fig.savefig(pdf_path)
        
        csv_path = pdf_path.replace(".pdf", ".csv")
        self.csv_fdroop(csv_path)
        
    def csv_fdroop(self, csv_path: Union[Path, str]):
        self.fdroop_log.to_csv(csv_path)
        
    def set_fdroop_log(self, df: pd.DataFrame):
        self.fdroop_log = df
        
    def __init__(self):
        self.fdroop_log = pd.DataFrame()
    
