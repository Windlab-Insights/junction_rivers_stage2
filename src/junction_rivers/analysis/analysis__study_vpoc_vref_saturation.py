import pandas as pd
from collections import OrderedDict
import signal_analysis
from standard_mains import standard_analysis_main, find_filepath_beneath_dir
from project_common import check_all_tx_settled, check_wtg_no_trip, check_bess_no_trip, check_wtg_no_frt, \
    check_bess_no_frt, check_vmp_no_frt
import json
import numpy as np

"""
Vref Steps Analysis
"""


def per_scenario_analysis(tuning: dict, scenario: dict, run_data: pd.DataFrame):
    init_vpoc_pu = float(scenario['Init_Vpoc_pu'])
    init_qpoc_pu = float(scenario['Init_Qpoc_pu'])
    init_ppoc_wind_mw = float(scenario['Init_Pwind_MW'])
    init_ppoc_bess_mw = float(scenario['Init_Pbess_MW'])
    init_vref_pu = float(scenario['Init_Vref_pu'])

    time_steps = json.loads(scenario['Time_Steps'])
    q_ctrl_mode = scenario['Control_Mode']

    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])
    post_init_time = model_init_time + float(scenario['Post_Init_Duration_s'])
    analysis_map = OrderedDict()

    qbase_mvar = float(tuning['SYS_Qbase_MVAr'])
    pbase_mw = float(tuning['SYS_Pbase_MW'])
    droop_ratio = float(tuning['SYS_Vdroop_perc'])

    # Constants that determine settling accuracy!
    min_settle_tol = (50 / 1000)
    vpoc_pu_tol = 0.01
    ppoc_tol_mw = 0.01 * pbase_mw
    qpoc_tol_mvar = 0.01 * qbase_mvar

    data = run_data[model_init_time:post_init_time]

    vref_step_start = time_steps[1] + model_init_time
    vref_step_end = time_steps[2] + model_init_time

    print(vref_step_start, vref_step_end)

    ppoc_mw = data["POC_P_MW"][model_init_time::]
    qpoc_mvar = data["POC_Q_MVAr"][model_init_time::]
    vpoc_pu = data["POC_V_pu"][model_init_time::]
    vref_pu = data["Vref_pu"][model_init_time::]
    if q_ctrl_mode == "VAR":
        qref_mvar = data["Qref_MVAr"][model_init_time::]
    else:
        qref_mvar = data["Qref_droop_MVAr"][model_init_time::]

    # droop_compensated_vref = np.clip((qpoc_mvar / qbase_mvar), -1, 1) * droop_ratio + vref_pu
    # vpoc_error_pu = droop_compensated_vref - vpoc_pu


    # Calculate all Settling Times pre- and post-step.
    pre_ppoc_settled, pre_ppoc_mw, _ = signal_analysis.calc_settled_value(ppoc_mw, vref_step_start, min_settle_tol,
                                                                          ppoc_tol_mw)
    post_ppoc_settled, post_ppoc_mw, _ = signal_analysis.calc_settled_value(ppoc_mw, vref_step_end, min_settle_tol,
                                                                            ppoc_tol_mw)

    pre_qpoc_settled, pre_qpoc_mvar, _ = signal_analysis.calc_settled_value(qpoc_mvar, vref_step_start, min_settle_tol,
                                                                            qpoc_tol_mvar)
    post_qpoc_settled, post_qpoc_mvar, _ = signal_analysis.calc_settled_value(qpoc_mvar, vref_step_end, min_settle_tol,
                                                                              qpoc_tol_mvar)

    pre_vpoc_settled, pre_vpoc_pu, _ = signal_analysis.calc_settled_value(vpoc_pu, vref_step_start, min_settle_tol,
                                                                          vpoc_pu_tol)
    post_vpoc_settled, post_vpoc_pu, post_vpoc_settled_time = signal_analysis.calc_settled_value(vpoc_pu, vref_step_end, min_settle_tol,
                                                                            vpoc_pu_tol)

    pre_qref_settled, pre_qref_mvar, _ = signal_analysis.calc_settled_value(qref_mvar, vref_step_start, min_settle_tol,
                                                                            4*qpoc_tol_mvar)
    post_qref_settled, post_qref_mvar, _ = signal_analysis.calc_settled_value(qref_mvar, vref_step_end, min_settle_tol,
                                                                              4*qpoc_tol_mvar)
    # Check if important signals have settled for this step.
    analysis_map[f"Ppoc_Settled"] = pre_ppoc_settled and post_ppoc_settled
    analysis_map[f"Qpoc_Settled"] = pre_qpoc_settled and post_qpoc_settled
    analysis_map[f"Vpoc_Settled"] = pre_vpoc_settled and post_vpoc_settled
    analysis_map[f"Qref_Settled"] = pre_qref_settled and post_qref_settled

    # Compute Rise/Settle time for Ppoc
    if analysis_map[f"Ppoc_Settled"]:
        analysis_map[f"Ppoc_Settled_MW"] = post_ppoc_mw
        t_settle = signal_analysis.calc_aemo_step_settling_time(ppoc_mw, vref_step_start, pre_ppoc_mw, vref_step_end,
                                                                post_ppoc_mw)
        analysis_map[f"Ppoc_AEMO_Settling_Time"] = t_settle
        if max(ppoc_mw[vref_step_start:vref_step_end]) - min(ppoc_mw[vref_step_start:vref_step_end]) < ppoc_tol_mw:
            analysis_map[f"Ppoc_AEMO_Settling_Time"] = 0


    # Compute Rise/Settle time for Qpoc
    if analysis_map[f"Qpoc_Settled"]:
        analysis_map[f"Qpoc_Settled_MVAr"] = post_qpoc_mvar
        analysis_map[f"Qpoc_Settled_pu"] = post_qpoc_mvar / qbase_mvar
        t_rise = signal_analysis.calc_aemo_step_rise_time(qpoc_mvar, vref_step_start, pre_qpoc_mvar, vref_step_end, post_qpoc_mvar)
        t_settle = signal_analysis.calc_aemo_step_settling_time(qpoc_mvar, vref_step_start, pre_qpoc_mvar, vref_step_end, post_qpoc_mvar)
        analysis_map[f"Qpoc_AEMO_Rise_Time"] = t_rise
        analysis_map[f"Qpoc_AEMO_Settling_Time"] = t_settle

    # Compute Rise/Settle time for Vpoc
    if analysis_map[f"Vpoc_Settled"]:
        analysis_map[f"Vpoc_Settled_pu"] = post_vpoc_pu
        t_settle = signal_analysis.calc_aemo_step_settling_time(vpoc_pu, vref_step_start, pre_vpoc_pu, vref_step_end, post_vpoc_pu)
        _, _, _, _, vpoc_ht = signal_analysis.calc_transient_characteristics(vpoc_pu, vref_step_start, pre_vpoc_pu, vref_step_end, post_vpoc_pu)
        analysis_map[f"Vpoc_AEMO_Settling_Time"] = t_settle
        analysis_map[f"Vpoc_Halving_Time"] = vpoc_ht

        if analysis_map[f"Qpoc_Settled"]:
            if abs(analysis_map[f"Qpoc_Settled_MVAr"] / qbase_mvar) < 0.99:
                vref_compensated = vref_pu - np.clip((qpoc_mvar / qbase_mvar), -1, 1) * droop_ratio
                verr = vpoc_pu - vref_compensated
                analysis_map[f"Vpoc_VrefDroop_pu"] = verr[t_settle::].values[0]
            else:
                analysis_map[f"Vpoc_VrefDroop_pu"] = -1

    analysis_map[f"Vref_Settled_pu"] = vref_pu[vref_step_start::].values[0]

    # Compute Rise/Settle time for Qpoc
    if analysis_map[f"Qref_Settled"] and analysis_map[f"Qpoc_Settled"]:
        analysis_map[f"Qpoc_to_Qref_Error_pu"] = abs(post_qref_mvar - post_qpoc_mvar) / qbase_mvar

    wt1_v_pu = data["WT1_V_pu"][vref_step_start:vref_step_end]
    wt2_v_pu = data["WT2_V_pu"][vref_step_start:vref_step_end]
    wt3_v_pu = data["WT3_V_pu"][vref_step_start:vref_step_end]
    wt4_v_pu = data["WT4_V_pu"][vref_step_start:vref_step_end]
    wt_worstcase_max = max(max(wt1_v_pu), max(wt2_v_pu), max(wt3_v_pu), max(wt4_v_pu))
    wt_worstcase_min = min(min(wt1_v_pu), min(wt2_v_pu), min(wt3_v_pu), min(wt4_v_pu))
    analysis_map[f"Max_WTG_V_pu"] = wt_worstcase_max
    analysis_map[f"Min_WTG_V_pu"] = wt_worstcase_min

    analysis_map['WTG_No_Trip'] = check_wtg_no_trip(data)
    analysis_map['BESS_No_Trip'] = check_bess_no_trip(data)
    analysis_map['VMP_No_FRT'] = check_vmp_no_frt(data)
    analysis_map['WTG_No_FRT'] = check_wtg_no_frt(data)
    analysis_map['BESS_No_FRT'] = check_bess_no_frt(data)

    return analysis_map


# def summary_analysis(tuning, analysis_df, data_root_dir, output_dir_path):


if __name__ == '__main__':
    analysis_df = standard_analysis_main(
        file_name_prefix=["Study_VpocSat", "Study_VrefSat"], per_scenario_analysis=per_scenario_analysis,
        summary_analysis=None)
