

import pandas as pd


def s52514_pref_step_analysis(spec_dict: dict, df: pd.DataFrame) -> dict : 
    
    analysis_map = {}
    
    nemde_interval_sec = float(spec_dict["NEMDE_Interval_sec_v"])
    init_time = float(spec_dict["substitutions"]["TIME_Full_Init_Time_sec"])
    
    poc_pref_mw = df["POC_WTG_Pref_MW"] # + df["POC_BESS_Pref_MW"]
    poc_p_mw = df["POC_P_MW"]
    
    previous_pref = poc_pref_mw.values[0]
    expected_poc_p_ramp_t = [0]
    expected_poc_p_ramp_mw = [previous_pref]
    
    for i, pref in zip(poc_pref_mw.index, poc_pref_mw.values):
        
        if pref != previous_pref:
            
            expected_poc_p_ramp_t.append(i-init_time)
            expected_poc_p_ramp_mw.append(previous_pref)
            
            expected_poc_p_ramp_t.append(i + nemde_interval_sec-init_time)
            expected_poc_p_ramp_mw.append(pref)
        
        previous_pref = pref
        
    expected_poc_p_ramp_t.append(i-init_time)
    expected_poc_p_ramp_mw.append(expected_poc_p_ramp_mw[-1])
    
    analysis_map = {
            "expected_poc_p_ramp_t" : expected_poc_p_ramp_t,
            "expected_poc_p_ramp_mw" : expected_poc_p_ramp_mw,
        }
    
    return analysis_map