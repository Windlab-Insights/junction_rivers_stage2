import pandas as pd
from collections import OrderedDict
from  bungaban.analysis import signal_analysis
from  bungaban.analysis.standard_mains import standard_analysis_main, find_filepath_beneath_dir
from  bungaban.analysis.project_common import check_all_tx_settled, check_wtg_no_trip, check_bess_no_trip, check_wtg_no_hvrt, check_bess_no_frt, check_vmp_no_frt
import numpy as np


"""
High-Voltage Faults.
"""
#

def per_scenario_analysis(tuning: dict, scenario: dict, run_data: pd.DataFrame):
    init_vpoc_pu = float(scenario['Init_Vpoc_pu'])
    init_qpoc_pu = float(scenario['Init_Qpoc_pu'])
    init_ppoc_wind_mw = float(scenario['Init_Pwind_MW'])
    init_ppoc_bess_mw = float(scenario['Init_Pbess_MW'])
    init_vref_pu = float(scenario['Init_Vref_pu'])

    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])
    post_init_time = model_init_time + float(scenario['Post_Init_Duration_s'])

    fault_start_time = model_init_time + float(scenario['Fault_Time'])
    fault_end_time = fault_start_time + float(scenario['Fault_Duration'])
    post_fault_end_time = fault_end_time + float(scenario['Post_Fault_Duration'])

    data = run_data[model_init_time:post_init_time]
    analysis_map = OrderedDict()

    # Constants that determine settling accuracy!
    min_settle_tol = (50 / 1000)
    ipoc_pu_tol = 0.07
    vpoc_pu_tol = 0.07

    analysis_map['WTG_No_Trip'] = check_wtg_no_trip(data)
    analysis_map['BESS_No_Trip'] = check_bess_no_trip(data)
    # analysis_map['VMP_No_FRT'] = check_vmp_no_frt(data)
    analysis_map['WTG_No_FRT'] = check_wtg_no_hvrt(data)
    # analysis_map['BESS_No_FRT'] = check_bess_no_frt(data)

    # Not a valid test if a trip occurs.
    valid_test = analysis_map['WTG_No_Trip'] and analysis_map['BESS_No_Trip']

    # Compute WTG Iq Values.
    all_wtg_settled = True
    for i in [1, 2, 3, 4, 5, 6]:
        iq_pos_pu = data[f'WTG_Iq_pu:1:{i}']
        v_pos_pu = data[f'WTG_Terminal_V_pu:1:{i}']
        pre_fault_iq_settled, pre_fault_iq_settled_value, _ = \
            signal_analysis.calc_settled_value(iq_pos_pu, fault_start_time, min_settle_tol, ipoc_pu_tol)
        fault_iq_settled, fault_iq_settled_value, _ = \
            signal_analysis.calc_settled_value(iq_pos_pu, fault_end_time, min_settle_tol, ipoc_pu_tol)
        fault_v_settled, fault_v_settled_value, _ = \
            signal_analysis.calc_settled_value(v_pos_pu, fault_end_time, min_settle_tol, vpoc_pu_tol)

        analysis_map[f"WT{i}_Signals_Settled"] = pre_fault_iq_settled and fault_iq_settled and fault_v_settled
        if analysis_map[f"WT{i}_Signals_Settled"]:
            analysis_map[f"WT{i}_diq_pu"] = fault_iq_settled_value - pre_fault_iq_settled_value
            analysis_map[f"WT{i}_Fault_iq_pu"] = fault_iq_settled_value
            analysis_map[f"WT{i}_Fault_V_pu"] = fault_v_settled_value

        all_wtg_settled = all_wtg_settled and analysis_map[f"WT{i}_Signals_Settled"]


    # Compute Bess Iq Values.
    all_bess_settled = True
    for i in [1, 2]:
        iq_pos_pu = data[f'BESS_Iq_pu:1:{i}']
        v_pos_pu = data[f'BESS_Terminal_V_pu:1:{i}']
        pre_fault_iq_settled, pre_fault_iq_settled_value, _ = \
            signal_analysis.calc_settled_value(iq_pos_pu, fault_start_time, min_settle_tol, ipoc_pu_tol)
        fault_iq_settled, fault_iq_settled_value, _ = \
            signal_analysis.calc_settled_value(iq_pos_pu, fault_end_time, min_settle_tol, ipoc_pu_tol)
        fault_v_settled, fault_v_settled_value, _ = \
            signal_analysis.calc_settled_value(v_pos_pu, fault_end_time, min_settle_tol, vpoc_pu_tol)

        analysis_map[f"BESS{i}_Signals_Settled"] = pre_fault_iq_settled and fault_iq_settled and fault_v_settled
        if analysis_map[f"BESS{i}_Signals_Settled"]:
            analysis_map[f"BESS{i}_diq_pu"] = fault_iq_settled_value - pre_fault_iq_settled_value
            analysis_map[f"BESS{i}_Fault_iq_pu"] = fault_iq_settled_value
            analysis_map[f"BESS{i}_Fault_V_pu"] = fault_v_settled_value
        all_bess_settled = all_bess_settled and analysis_map[f"BESS{i}_Signals_Settled"]


    # Compute POC Iq Values.
    poc_iq_pos_pu = data['POC_Iq_pos_pu']
    poc_i_pu = data['POC_I_kA'] / float(tuning['SYS_Ibase_kA'])

    fault_ipoc_settled, fault_ipoc_settled_value, _ = \
        signal_analysis.calc_settled_value(poc_i_pu, fault_end_time, min_settle_tol, ipoc_pu_tol)
    valid_test = valid_test and fault_ipoc_settled
    analysis_map["Fault_Ipoc_mcc_pu"] = fault_ipoc_settled_value

    pre_fault_iq_settled, pre_fault_iq_settled_value, _ = \
        signal_analysis.calc_settled_value(poc_iq_pos_pu, fault_start_time, min_settle_tol, ipoc_pu_tol)
    fault_iq_settled, fault_iq_settled_value, _ = \
        signal_analysis.calc_settled_value(poc_iq_pos_pu, fault_end_time, min_settle_tol, ipoc_pu_tol)



    valid_test = valid_test and pre_fault_iq_settled and fault_iq_settled and all_wtg_settled and all_bess_settled

    analysis_map["Prefault_iq_settled"] = pre_fault_iq_settled
    analysis_map["Prefault_iq_pu"] = pre_fault_iq_settled_value

    analysis_map["Fault_iq_settled"] = fault_iq_settled
    analysis_map["Fault_iq_pu"] = fault_iq_settled_value
    analysis_map["Fault_iq_mcc_pu"] = fault_iq_settled_value
    analysis_map["Delta_iq_mcc_pu"] = (fault_iq_settled_value - pre_fault_iq_settled_value)

    # Compute AEMO Rise and Settling time POC Iq if the signal was settled.
    if pre_fault_iq_settled and fault_iq_settled:
        input_pars = (poc_iq_pos_pu, fault_start_time, pre_fault_iq_settled_value, fault_end_time,
                      fault_iq_settled_value)
        iq_rise_time = signal_analysis.calc_aemo_step_rise_time(*input_pars)
        iq_settle_time = signal_analysis.calc_aemo_step_settling_time(*input_pars)
        mas_adeq_damped_diff, aas_adeq_damp_diff = signal_analysis.calc_aemo_adequately_damped(*input_pars)

        analysis_map["Fault_iq_AEMO_Rise_Time"] = iq_rise_time
        analysis_map["Fault_iq_AEMO_Settling_Time"] = iq_settle_time
        analysis_map["Fault_iq_AEMO_MAS_Adeq_Damp"] = mas_adeq_damped_diff >= 0
        analysis_map["Fault_iq_AEMO_AAS_Adeq_Damp"] = aas_adeq_damp_diff >= 0

    vpoc_pu = data["POC_V_pu"]

   # Compute POC Voltage Values.
    pre_fault_vpoc_settled, pre_fault_vpoc_settled_value, _ = \
        signal_analysis.calc_settled_value(vpoc_pu, fault_start_time, min_settle_tol, vpoc_pu_tol)
    fault_vpoc_settled, fault_vpoc_settled_value, _ = \
        signal_analysis.calc_settled_value(vpoc_pu, fault_end_time, min_settle_tol, vpoc_pu_tol)

    valid_test = valid_test and pre_fault_vpoc_settled and fault_vpoc_settled
    analysis_map["Prefault_Vpoc_settled"] = pre_fault_vpoc_settled
    analysis_map["Prefault_Vpoc_pu"] = pre_fault_vpoc_settled_value
    analysis_map["Fault_Vpoc_settled"] = fault_vpoc_settled
    analysis_map["Fault_Vpoc_pu"] = fault_vpoc_settled_value
    analysis_map["Delta_V_pu(rel. 120%)"] = fault_vpoc_settled_value - 1.20
    analysis_map["diq_mcc/dV"] = analysis_map["Delta_iq_mcc_pu"] / analysis_map["Delta_V_pu(rel. 120%)"]

    analysis_map["Fault_Vpoc_above_115"]= 1.15 < min(vpoc_pu[fault_start_time+0.05:fault_end_time])
    analysis_map["Fault_Vpoc_above_120"] = 1.20 < min(vpoc_pu[fault_start_time + 0.05:fault_end_time])

    # diff = vpoc_pu - data['WT_max_Vrms_All_Phases_pu']
    #
    # analysis_map["POC_WTG_min_Diff"] = min(diff)
    # analysis_map["min_WT_max_Vrms_All_Phases_pu"] = min(data['WT_max_Vrms_All_Phases_pu'])


    # Compute Post Fault Active-Power Recovery Time.
    # Find the time it takes for the signal POC_P_MW to remain above 0.95*initial_ppoc_mw

    # Init. Ppoc.
    initial_ppoc_mw = data["POC_P_MW"].values[0]
    sig = data["POC_P_MW"][fault_end_time:post_fault_end_time]
    if initial_ppoc_mw > 0:
        within_95pc_check = sig > 0.95 * initial_ppoc_mw
    else:
        within_95pc_check = sig < 0.95 * initial_ppoc_mw

    if list(sig[within_95pc_check].index) == []:
        analysis_map["Ppoc_First_95pc_Recovery_Time"] = None
    else:
        analysis_map["Ppoc_First_95pc_Recovery_Time"] = min(list(sig[within_95pc_check].index)) - fault_end_time

    within_95pc_check = within_95pc_check[::-1].iteritems()
    within_95pc_time = post_fault_end_time
    for curr_time, currently_above in within_95pc_check:
        if currently_above:
            within_95pc_time = curr_time
        else:
            break
    analysis_map["Ppoc_95pc_Recovery_Time"] = within_95pc_time - fault_end_time
    # ----

    # PPOC Within 10 MW
    if initial_ppoc_mw > 0:
        within_10mw_check = (sig > (initial_ppoc_mw - 10))
    else:
        within_10mw_check = (sig < (initial_ppoc_mw + 10))

    if list(sig[within_10mw_check].index) == []:
        analysis_map[f"Ppoc_First_Within10MW_Recovery_Time"] = None
    else:
        # print(list(sig[within_10mw_check].index))
        analysis_map[f"Ppoc_First_Within10MW_Recovery_Time"] = \
            min(list(sig[within_10mw_check].index)) - fault_end_time
    within_10mw_check = within_10mw_check[::-1].iteritems()
    within_10mw_time = post_fault_end_time
    for curr_time, currently_above in within_10mw_check:
        if currently_above:
            within_10mw_time = curr_time
        else:
            break
    analysis_map[f"Ppoc_Within10MW_Recovery_Time"] = within_10mw_time - fault_end_time

    # Iq AAS and MAS Settle-time. Not real settle time. Time till we are 'above' the mas/aas standard.
    mas_injection = np.clip(2 * (1.2 - fault_vpoc_settled_value), -1, 0)
    aas_injection = np.clip(6 * (1.15 - fault_vpoc_settled_value), -1, 0)
    max_cont_curr = 1
    mas_iq_target = np.clip(mas_injection + pre_fault_iq_settled_value, -max_cont_curr, max_cont_curr)
    aas_iq_target = np.clip(aas_injection + pre_fault_iq_settled_value, -max_cont_curr, max_cont_curr)
    sig = poc_iq_pos_pu[fault_start_time:fault_end_time]

    below_mas_check = (sig < mas_iq_target)[::-1].iteritems()
    below_mas_time = fault_end_time
    for curr_time, currently_above in below_mas_check:
        if currently_above:
            below_mas_time = curr_time
        else:
            break
    analysis_map["Iq_Below_MAS_Time"] = below_mas_time - fault_start_time

    below_aas_check = (sig < aas_iq_target)[::-1].iteritems()
    below_aas_time = fault_end_time
    for curr_time, currently_above in below_aas_check:
        if currently_above:
            below_aas_time = curr_time
        else:
            break
    analysis_map["Iq_Below_AAS_Time"] = below_aas_time - fault_start_time

    # Check if the windfarm tripped.
    analysis_map["Valid_Test"] = valid_test
    #
    # # Look over all FRT engages. and calculate lowest voltage.
    # for frt_sig in ['WT1_FRT', 'WT2_FRT', 'WT3_FRT', 'WT4_FRT']:
    #     hvrt_only = data[frt_sig].copy(deep=True)
    #     hvrt_only[data[frt_sig] == -1] = 0
    #
    #     if frt_sig == 'WT1_FRT':
    #         hvrt_on_filter = (hvrt_only.diff() == 1)
    #         hvrt_off_filter = (hvrt_only.diff() == -1)
    #         if len(hvrt_on_filter) > 0:
    #             analysis_map["Max. WT_V_pu(HVRT_on)"] = max(data[frt_sig][hvrt_on_filter].values)
    #         if len(hvrt_off_filter) > 0:
    #             analysis_map["Max. WT_V_pu(HVRT_off)"] = max(data[frt_sig][hvrt_off_filter].values)
    #     else:
    #         hvrt_on_filter = hvrt_on_filter | (hvrt_only.diff() == 1)
    #         hvrt_off_filter = hvrt_off_filter | (hvrt_only.diff() == -1)
    #         if len(hvrt_on_filter) > 0:
    #             analysis_map["Min. WT_V_pu(HVRT_on)"] = min(analysis_map["Min. WT_V_pu(HVRT_on)"],
    #                                                         min(data[frt_sig][hvrt_on_filter].values))
    #
    #
    # if len(data['POC_V_pu'][hvrt_on_filter]) > 0:
    #     analysis_map["Min. WT_V_pu(LVRT_on)"] = min(data['WT1_V_pu'][hvrt_on_filter].values)
    #     analysis_map["Min. POC_V_PU(LVRT_on)"] = min(data['POC_V_pu'][hvrt_on_filter].values)
    # if len(data['POC_V_pu'][hvrt_off_filter]) > 0:
    #     analysis_map["Min. WT_V_pu(LVRT_off)"] = min(data['WT1_V_pu'][hvrt_on_filter].values)
    #     analysis_map["Min. POC_V_PU(LVRT_off)"] = min(data['POC_V_pu'][hvrt_off_filter].values)

    return analysis_map


# def summary_analysis(tuning, analysis_df, data_root_dir, output_dir_path):


if __name__ == '__main__':
    analysis_df = standard_analysis_main(
        file_name_prefix=["s5255_Fault_TOV"], per_scenario_analysis=per_scenario_analysis,
        summary_analysis=None)

