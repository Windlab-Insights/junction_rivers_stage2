import pandas as pd
from icecream import ic
from collections import OrderedDict
import signal_analysis
from standard_mains import standard_analysis_main

"""
Freqency Disturbance Analysis
"""

def per_scenario_analysis(tuning, scenario: pd.DataFrame, run_data: pd.DataFrame):


    init_vpoc_pu = float(scenario['Init_Vpoc_pu'])
    init_qpoc_pu = float(scenario['Init_Qpoc_pu'])
    init_ppoc_wind_mw = float(scenario['Pref_Wind_MW_v'])
    init_ppoc_bess_mw = float(scenario['Pref_BESS_MW_v'])
    init_vref_pu = float(scenario['Vref_pu_v'])

    ic(init_vpoc_pu, init_qpoc_pu, init_ppoc_wind_mw, init_ppoc_bess_mw)

    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])
    post_init_time = model_init_time + float(scenario['Post_Init_Duration_s'])

    data = run_data[model_init_time:post_init_time]

    ppoc_mw = data["POC_P_MW"]
    qpoc_mvar = data["POC_Q_MVAr"]
    vpoc_pu = data["POC_V_pu"]
    vref_pu = data["Vref_pu"]
    poc_freq_hz =  data["Fpoc_Hz"]

    ppoc_target_mw = init_ppoc_wind_mw + init_ppoc_bess_mw

    analysis_map = OrderedDict()

    start_of_dist_s = 1.0  # scenario["Grid_Hz_t"].split(",")[1]
    end_of_dist_t = 13.0
    sample_time = 12.0

    freq_during_disturbance_hz = poc_freq_hz.values[-1] #[start_of_dist_s:end_of_dist_t].mean()
    poc_p_during_disturbance_mw = ppoc_mw.values[-1] #[start_of_dist_s:end_of_dist_t].mean()
    wtg_p_during_disturbance_mw = 15*data["WT1_P_MW"].values[-1] + 17*data["WT2_P_MW"].values[-1] + 18*data["WT3_P_MW"].values[-1] + 18*data["WT4_P_MW"].values[-1]
    bess_p_during_disturbance_mw = data["BESS1_P_MW"].values[-1] + data["BESS2_P_MW"].values[-1] + data["BESS3_P_MW"].values[-1] + data["BESS4_P_MW"].values[-1]




    analysis_map["Freq_During_Dist_Hz"] = freq_during_disturbance_hz
    analysis_map["POC_P_During_Dist_MW"] = poc_p_during_disturbance_mw
    analysis_map["WTG_P_During_Dist_MW"] = wtg_p_during_disturbance_mw
    analysis_map["BESS_P_During_Dist_MW"] = bess_p_during_disturbance_mw

    analysis_map["Delta_Freq_During_Dist_Hz"] = freq_during_disturbance_hz - 50
    analysis_map["Delta_POC_P_During_Dist_MW"] = poc_p_during_disturbance_mw - ppoc_target_mw
    analysis_map["Delta_WTG_P_During_Dist_MW"] = wtg_p_during_disturbance_mw - init_ppoc_wind_mw
    analysis_map["Delta_BESS_P_During_Dist_MW"] = bess_p_during_disturbance_mw - init_ppoc_wind_mw

    return analysis_map

#
#     analysis_map['Ppoc_Variation'] = abs(ppoc_mw.max() - ppoc_mw.min())
#     analysis_map['Ppoc_Average'] = ppoc_mw.mean()
#     analysis_map['Ppoc_Accurate'] = signal_analysis.signal_is_constant(raw_signal=ppoc_mw,
#                                        value=ppoc_target_mw,
#                                        variation_tolerance=0.02 * float(tuning['SYS_Pbase_MW']),
#                                        value_tolerance=0.02 * float(tuning['SYS_Pbase_MW']))
#
#     analysis_map['Vpoc_Variation'] = abs(vpoc_pu.max() - vpoc_pu.min())
#     analysis_map['Vpoc_Average'] = vpoc_pu.mean()
#     analysis_map['Vpoc_Accurate'] = signal_analysis.signal_is_constant(raw_signal=vpoc_pu,
#                                        value=init_vpoc_pu,
#                                        variation_tolerance=0.01,
#                                        value_tolerance=0.01)
#
#     if scenario['Control_Mode'].lower() == 'droop':
#         import numpy as np
#         droop_qref_pu = np.clip((init_vref_pu - init_vpoc_pu)/float(tuning['SYS_Vdroop_perc']), -1, 1)
#         qbase_mvar = float(tuning['SYS_Qbase_MVAr'])
#         droop_qref_mvar = droop_qref_pu * qbase_mvar
#         analysis_map['Qpoc_Variation'] = abs(qpoc_mvar.max() - qpoc_mvar.min())
#         analysis_map['Qpoc_Average'] = qpoc_mvar.mean()
#         analysis_map['Qpoc_Accurate'] = signal_analysis.signal_is_constant(raw_signal=qpoc_mvar,
#                                            value=droop_qref_mvar,
#                                            variation_tolerance=0.01 * qbase_mvar,
#                                            value_tolerance=0.01 * qbase_mvar)
#
#     tx_settled = []
#     wtg_no_trip = []
#     wtg_no_frt = []
#     bess_no_trip = []
#     bess_no_frt = []
#
#     for idx in range(1, 5):
#         tx_settled.append(signal_analysis.signal_is_constant(data[f'TX{idx}_Settled'], value=1))
#         wtg_no_trip.append(signal_analysis.signal_is_constant(data[f'WT{idx}_Trip'], value=0))
#         bess_no_trip.append(signal_analysis.signal_is_constant(data[f'BESS{idx}_Trip'], value=0))
#         wtg_no_frt.append(signal_analysis.signal_is_constant(data[f'WT{idx}_FRT'], value=0))
#
#         bess_no_frt.append(signal_analysis.signal_is_constant(data[f'BESS{idx}_AVR_FRZ'], value=0))
#         bess_no_frt.append(signal_analysis.signal_is_constant(data[f'BESS{idx}_GOV_FRZ'], value=0))
#
#     analysis_map['TX_Settled'] = all(tx_settled)
#     analysis_map['WTG_No_Trip'] = all(wtg_no_trip)
#     analysis_map['BESS_No_Trip'] = all(bess_no_trip)
#     analysis_map['VMP_No_FRT'] = signal_analysis.signal_is_constant(data[f'VMP_FRT'], value=0)
#     analysis_map['WTG_No_FRT'] = all(wtg_no_frt)
#     analysis_map['BESS_No_FRT'] = all(bess_no_frt)
#
#     analysis_map['Compliance'] = False
#     if all(analysis_map[cond] for cond in ['Ppoc_Accurate', 'Qpoc_Accurate', 'Vpoc_Accurate',
#          'TX_Settled', 'WTG_No_Trip', 'BESS_No_Trip',
#          'VMP_No_FRT', 'WTG_No_FRT', 'BESS_No_FRT']):
#         analysis_map['Compliance'] = True
#



if __name__ == '__main__':
    analysis_df = standard_analysis_main(
        file_name_prefix=["s52511_GridHz"], 
        per_scenario_analysis=per_scenario_analysis,
        summary_analysis=None
    )

