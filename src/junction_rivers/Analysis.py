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
        analysis_map['Delta_BESS_P_During_Dist_MW'] = (bess_p_during_disturbance_mw - init_ppoc_bess_mw)/150
        
        self.fdroop_log = self.fdroop_log.append(other=pd.DataFrame([analysis_map]), ignore_index=True)

        return

    def plot_fdroop(self, pdf_path: Union[Path, str]):
        plt.clf()
        fig, axs = plt.subplots(2,2)
        ax1 = axs[0,0]
        ax2 = axs[0,1]
        ax3 = axs[1,0]
        ax4 = axs[1,1]
        ax1.set_title('Droop response large scale')
        ax1.set_xlabel('Frequency (Hz)')
        ax1.set_ylabel('Active Power (MW)')
        
        ax2.set_title('Droop response in per unit large scale')
        ax2.set_xlabel('Delta Frequency (pu)')
        ax2.set_ylabel('Delta Active Power (pu)')
        
        ax3.set_title('Droop response small scale')
        ax3.set_xlabel('Frequency (Hz)')
        ax3.set_ylabel('Active Power (MW)')
        
        ax4.set_title('Droop response in per unit small scale')
        ax4.set_xlabel('Delta Frequency (pu)')
        ax4.set_ylabel('Delta Active Power (pu)')
        
        
        fdroop_log_large = self.fdroop_log[(self.fdroop_log['Delta_Freq_During_Dist_Hz'] >= 0.004) | (self.fdroop_log['Delta_Freq_During_Dist_Hz'] <= -0.004)]
        fdroop_log_small = self.fdroop_log[(self.fdroop_log['Delta_Freq_During_Dist_Hz'] <= 0.004) & (self.fdroop_log['Delta_Freq_During_Dist_Hz'] >= -0.004)]
        
        l_freq_during_disturbance_hz = fdroop_log_large['Freq_During_Dist_Hz']
        l_poc_p_during_disturbance_mw = fdroop_log_large['POC_P_During_Dist_MW']
        l_wtg_p_during_disturbance_mw = fdroop_log_large['WTG_P_During_Dist_MW']
        l_bess_p_during_disturbance_mw = fdroop_log_large['BESS_P_During_Dist_MW']
        
        l_delta_freq_during_disturbance_hz = fdroop_log_large['Delta_Freq_During_Dist_Hz']
        l_delta_poc_p_during_disturbance_mw = fdroop_log_large['Delta_POC_P_During_Dist_MW']
        l_delta_wtg_p_during_disturbance_mw = fdroop_log_large['Delta_WTG_P_During_Dist_MW']
        l_delta_bess_p_during_disturbance_mw = fdroop_log_large['Delta_BESS_P_During_Dist_MW']
        
        s_freq_during_disturbance_hz = fdroop_log_small['Freq_During_Dist_Hz']
        s_poc_p_during_disturbance_mw = fdroop_log_small['POC_P_During_Dist_MW']
        s_wtg_p_during_disturbance_mw = fdroop_log_small['WTG_P_During_Dist_MW']
        s_bess_p_during_disturbance_mw = fdroop_log_small['BESS_P_During_Dist_MW']
        
        s_delta_freq_during_disturbance_hz = fdroop_log_small['Delta_Freq_During_Dist_Hz']
        s_delta_poc_p_during_disturbance_mw = fdroop_log_small['Delta_POC_P_During_Dist_MW']
        s_delta_wtg_p_during_disturbance_mw = fdroop_log_small['Delta_WTG_P_During_Dist_MW']
        s_delta_bess_p_during_disturbance_mw = fdroop_log_small['Delta_BESS_P_During_Dist_MW']
        
        ax1.plot(l_freq_during_disturbance_hz, l_poc_p_during_disturbance_mw, 'r', ls='', marker='o', label='POC')
        ax1.plot(l_freq_during_disturbance_hz, l_wtg_p_during_disturbance_mw, 'b', ls='', marker='o', label='WT')
        ax1.plot(l_freq_during_disturbance_hz, l_bess_p_during_disturbance_mw, 'g', ls='', marker='o', label='BESS')
        ax1.legend(loc="upper left")
        
        ax2.plot(l_delta_freq_during_disturbance_hz, l_delta_poc_p_during_disturbance_mw,'r', ls='', marker='o', label='POC')
        ax2.plot(l_delta_freq_during_disturbance_hz, l_delta_wtg_p_during_disturbance_mw, 'b', ls='', marker='o', label='WT')
        ax2.plot(l_delta_freq_during_disturbance_hz, l_delta_bess_p_during_disturbance_mw, 'g', ls='', marker='o', label='BESS')
        
        ax3.plot(s_freq_during_disturbance_hz, s_poc_p_during_disturbance_mw, 'r', ls='', marker='o', label='POC')
        ax3.plot(s_freq_during_disturbance_hz, s_wtg_p_during_disturbance_mw, 'b', ls='', marker='o', label='WT')
        ax3.plot(s_freq_during_disturbance_hz, s_bess_p_during_disturbance_mw, 'g', ls='', marker='o', label='BESS')
        
        ax4.plot(s_delta_freq_during_disturbance_hz, s_delta_poc_p_during_disturbance_mw,'r', ls='', marker='o', label='POC')
        ax4.plot(s_delta_freq_during_disturbance_hz, s_delta_wtg_p_during_disturbance_mw, 'b', ls='', marker='o', label='WT')
        ax4.plot(s_delta_freq_during_disturbance_hz, s_delta_bess_p_during_disturbance_mw, 'g', ls='', marker='o', label='BESS')
        
        
        fig.savefig(pdf_path)
        
        csv_path = pdf_path.replace(".pdf", ".csv")
        self.csv_fdroop(csv_path)
        
    def csv_fdroop(self, csv_path: Union[Path, str]):
        self.fdroop_log.to_csv(csv_path)
        
    def set_fdroop_log(self, df: pd.DataFrame):
        self.fdroop_log = df
        
    def __init__(self):
        self.fdroop_log = pd.DataFrame()
    
