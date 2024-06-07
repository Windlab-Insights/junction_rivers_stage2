import pandas as pd
from icecream import ic
from collections import OrderedDict
import signal_analysis
from standard_mains import standard_analysis_main
from project_common import qref_mvar_via_voltage_droop, check_all_tx_settled, check_wtg_no_trip, check_bess_no_trip, check_wtg_no_frt, \
    check_bess_no_frt, check_vmp_no_frt

"""
Flatrun Analysis
"""
#

def per_scenario_analysis(tuning: dict, scenario: dict, run_data: pd.DataFrame):
    init_vpoc_pu = float(scenario['Init_Vpoc_pu'])
    init_qpoc_pu = float(scenario['Init_Qpoc_pu'])
    init_ppoc_wind_mw = float(scenario['Init_Pwind_MW'])
    init_ppoc_bess_mw = float(scenario['Init_Pbess_MW'])
    init_vref_pu = float(scenario['Init_Vref_pu'])

    ic(init_vpoc_pu, init_qpoc_pu, init_ppoc_wind_mw, init_ppoc_bess_mw)

    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])
    post_init_time = model_init_time + float(scenario['Post_Init_Duration_s'])

    qbase_mvar = float(tuning['SYS_Qbase_MVAr'])

    data = run_data[model_init_time:post_init_time]

    ppoc_mw = data["POC_P_MW"]
    qpoc_mvar = data["POC_Q_MVAr"]
    vpoc_pu = data["POC_V_pu"]
    vref_pu = data["Vref_pu"]

    ppoc_target_mw = init_ppoc_wind_mw + init_ppoc_bess_mw
    analysis_map = OrderedDict()

    analysis_map['Ppoc_Variation'] = abs(ppoc_mw.max() - ppoc_mw.min())
    analysis_map['Ppoc_Average'] = ppoc_mw.mean()
    analysis_map['Ppoc_Accurate'] = signal_analysis.signal_is_constant(raw_signal=ppoc_mw,
                                       value=ppoc_target_mw,
                                       variation_tolerance=0.02 * float(tuning['SYS_Pbase_MW']),
                                       value_tolerance=0.02 * float(tuning['SYS_Pbase_MW']))

    analysis_map['Vpoc_Variation'] = abs(vpoc_pu.max() - vpoc_pu.min())
    analysis_map['Vpoc_Average'] = vpoc_pu.mean()
    analysis_map['Vpoc_Accurate'] = signal_analysis.signal_is_constant(raw_signal=vpoc_pu,
                                       value=init_vpoc_pu,
                                       variation_tolerance=0.01,
                                       value_tolerance=0.01)

    if scenario['Control_Mode'].lower() == 'droop':
        import numpy as np
        droop_qref_mvar = qref_mvar_via_voltage_droop(tuning, init_vref_pu, init_vpoc_pu)
        analysis_map['Qpoc_Variation'] = abs(qpoc_mvar.max() - qpoc_mvar.min())
        analysis_map['Qpoc_Average'] = qpoc_mvar.mean()
        analysis_map['Qpoc_Accurate'] = signal_analysis.signal_is_constant(raw_signal=qpoc_mvar,
                                           value=droop_qref_mvar,
                                           variation_tolerance=0.01 * qbase_mvar,
                                           value_tolerance=0.01 * qbase_mvar)


    analysis_map['TX_Settled'] = check_all_tx_settled(data)
    analysis_map['WTG_No_Trip'] = check_wtg_no_trip(data)
    analysis_map['BESS_No_Trip'] = check_bess_no_trip(data)
    analysis_map['VMP_No_FRT'] = check_vmp_no_frt(data)
    analysis_map['WTG_No_FRT'] = check_wtg_no_frt(data)
    analysis_map['BESS_No_FRT'] = check_bess_no_frt(data)

    analysis_map['Compliance'] = False
    if all(analysis_map[cond] for cond in ['Ppoc_Accurate', 'Qpoc_Accurate', 'Vpoc_Accurate',
         'TX_Settled', 'WTG_No_Trip', 'BESS_No_Trip',
         'VMP_No_FRT', 'WTG_No_FRT', 'BESS_No_FRT']):
        analysis_map['Compliance'] = True




    return analysis_map


# def summary_analysis(project_tuning, filtered_df, data_root_dir, output_dir_path):
# #     empty_dataframe = pd.DataFrame()
# #     print("in here")
# #     return empty_dataframe

if __name__ == '__main__':
    # analysis_df = standard_analysis_main(
    #     file_name_prefix=["s5251_FlatRun"], per_scenario_analysis=per_scenario_analysis,
    #     summary_analysis=summary_analysis)
    analysis_df = standard_analysis_main(
        file_name_prefix=["s5251_FlatRun"], per_scenario_analysis=per_scenario_analysis, summary_analysis=None)