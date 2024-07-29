

import pandas as pd


def s5253_s5258_freq_dist_analysis(spec_dict: dict, df: pd.DataFrame) -> dict : 
    
    over_freq_threshold = 52 # Replace with GS
    
    analysis_map = {}
    
    # Expected Over Frequency Active Power Ramp Down 
    if max(df["POC_Freq_Hz"].values) >= over_freq_threshold:
        
        active_power_rampdown_threshold_mw = (spec_dict["Init_Pwind_MW"] + spec_dict["Init_Pbess_MW"])/2
        
        poc_freq_hz = df["POC_Freq_Hz"]
        start_time_sec = float(spec_dict["substitutions"]["TIME_Full_Init_Time_sec"])
        active_power_rampdown_start_s = min(list(poc_freq_hz[poc_freq_hz > over_freq_threshold].index)) - start_time_sec
        
        active_power_rampdown_expected_end_s = active_power_rampdown_start_s + 7
        
        analysis_map = {
            "active_power_rampdown_threshold_mw" : active_power_rampdown_threshold_mw,
            "active_power_rampdown_start_s" : active_power_rampdown_start_s,
            "active_power_rampdown_expected_end_s" : active_power_rampdown_expected_end_s,
        }
        
    return analysis_map