import pandas as pd
from collections import OrderedDict
from  bungaban.analysis import signal_analysis
from bungaban.analysis.standard_mains import standard_analysis_main, find_filepath_beneath_dir
from bungaban.analysis.project_common import check_all_tx_settled, check_wtg_no_trip, check_bess_no_trip, check_wtg_no_lvrt, check_bess_no_frt, check_vmp_no_frt
import matplotlib.pyplot as plt

"""
Low-Voltage Faults.
"""

class PscadFaultType:
    NO_FAULT = 0
    A_TO_G = 1
    B_TO_G = 2
    C_TO_G = 3
    AB_TO_G = 4
    AC_TO_G = 5
    BC_TO_G = 6
    ABC_TO_G = 7
    A_TO_B = 8
    A_TO_C = 9
    B_TO_C = 10
    A_TO_B_TO_C = 11


def healthy_voltage_map(data):
    healthy_lookup = {PscadFaultType.ABC_TO_G: {"min_fault_pu": None, "max_healthy_pu": None},
                      PscadFaultType.A_TO_B_TO_C: {"min_fault_pu": None, "max_healthy_pu": None}}
    all_phases = set(["A", "B", "C"])

    for fault_type_no in [PscadFaultType.ABC_TO_G, PscadFaultType.A_TO_B_TO_C]:
        healthy_phases = list(all_phases)
        min_faulted_phase = min(min(data[f"POC_PhaseA_V_pu"]),
                                min(data[f"POC_PhaseB_V_pu"]),
                                min(data[f"POC_PhaseC_V_pu"]))

        max_healthy_phase = None
        # max_healthy_phase = max(max(data[f"POC_PhaseA_V_pu"]),
        #                         max(data[f"POC_PhaseB_V_pu"]),
        #                         max(data[f"POC_PhaseC_V_pu"]))
        healthy_lookup[fault_type_no] = {"min_fault_pu": min_faulted_phase, "max_healthy_pu": max_healthy_phase}

    for faulted_phase, fault_type_no in [("A", PscadFaultType.A_TO_G),
                                         ("B", PscadFaultType.B_TO_G),
                                         ("C", PscadFaultType.C_TO_G)]:
        healthy_phases = list(all_phases - set([faulted_phase]))
        min_faulted_phase = min(data[f"POC_Phase{faulted_phase}_V_pu"])
        max_healthy_phase = max(max(data[f"POC_Phase{healthy_phases[0]}_V_pu"]),
                                max(data[f"POC_Phase{healthy_phases[1]}_V_pu"]))
        healthy_lookup[fault_type_no] = {"min_fault_pu": min_faulted_phase, "max_healthy_pu": max_healthy_phase}

    for faulted_phases, llg_no, ll_no in [(("A", "B"), PscadFaultType.AB_TO_G, PscadFaultType.A_TO_B),
                                          (("B", "C"), PscadFaultType.BC_TO_G, PscadFaultType.B_TO_C),
                                          (("A", "C"), PscadFaultType.AC_TO_G, PscadFaultType.A_TO_C)]:
        healthy_phase = list(all_phases - set(faulted_phases))[0]
        min_faulted_phase = min(min(data[f"POC_Phase{faulted_phases[0]}_V_pu"]),
                                min(data[f"POC_Phase{faulted_phases[1]}_V_pu"]))
        max_healthy_phase = max(data[f"POC_Phase{healthy_phase}_V_pu"])
        result = {"min_fault_pu": min_faulted_phase, "max_healthy_pu": max_healthy_phase}
        healthy_lookup[llg_no] = result
        healthy_lookup[ll_no] = result
    return healthy_lookup


def per_scenario_analysis(tuning: dict, scenario: dict, run_data: pd.DataFrame):
    init_vpoc_pu = float(scenario['Init_Vpoc_pu'])
    init_qpoc_pu = float(scenario['Init_Qpoc_pu'])
    init_ppoc_wind_mw = float(scenario['Init_Pwind_MW'])
    init_ppoc_bess_mw = float(scenario['Init_Pbess_MW'])
    init_vref_pu = float(scenario['Init_Vref_pu'])
    fault_type = int(scenario['Fault_Type'])

    valid_test = True

    # S5.1a.8
    max_clearance_time_secs = 9999

    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])
    post_init_time = model_init_time + float(scenario['Post_Init_Duration_s'])

    fault_start_time = model_init_time + float(scenario['Fault_Time'])
    fault_end_time = fault_start_time + min(float(scenario['Fault_Duration']), max_clearance_time_secs)
    post_fault_start_time = fault_start_time + float(scenario['Fault_Duration'])
    post_fault_end_time = fault_end_time + float(scenario['Post_Fault_Duration'])

    data = run_data[model_init_time:post_init_time]
    analysis_map = OrderedDict()

    # Constants that determine settling accuracy!
    min_settle_tol = (50 / 1000)
    ipoc_pu_tol = 0.02
    vpoc_pu_tol = 0.02

    # # analysis_map['WTG_No_Trip'] = check_wtg_no_trip(data)
    # # analysis_map['BESS_No_Trip'] = check_bess_no_trip(data)
    # analysis_map['VMP_No_FRT'] = check_vmp_no_frt(data)
    # analysis_map['WTG_No_FRT'] = check_wtg_no_frt(data)
    # analysis_map['BESS_No_FRT'] = check_bess_no_frt(data)

    # Not a valid test if a trip occurs.
    # # valid_test = analysis_map['WTG_No_Trip'] and analysis_map['BESS_No_Trip']
    fault_data = data[fault_start_time:fault_end_time]
    below_80pc = fault_data[fault_data[f'POC_Vrms_pu'] < 0.8]

    if len(below_80pc) > 0:
        detected_fault_start_time = below_80pc.index[0]
        # print(fault_start_time, detected_fault_start_time)
    else:
        detected_fault_start_time = fault_start_time

    # Compute WTG Iq and P Values.
    wtg_max_iq_rise_time = 0
    wtg_max_iq_settle_time = 0
    all_wtg_settled = True
    for i in [1, 2, 3, 4, 5, 6]:
        iq_pos_pu = data[f'WTG_Iq_pu:{i}']
        v_pos_pu = data[f'WTG_MV_V_pu:{i}']
        pre_fault_iq_settled, pre_fault_iq_settled_value, _ = \
            signal_analysis.calc_settled_value(iq_pos_pu, fault_start_time, min_settle_tol, ipoc_pu_tol) 
        fault_iq_settled, fault_iq_settled_value, fault_iq_settled_time = \
            signal_analysis.calc_settled_value(iq_pos_pu, fault_end_time, min_settle_tol, ipoc_pu_tol)
        fault_v_settled, fault_v_settled_value, _ = \
            signal_analysis.calc_settled_value(v_pos_pu, fault_end_time, min_settle_tol, vpoc_pu_tol)

        analysis_map[f"WT{i}_Signals_Settled"] = pre_fault_iq_settled and fault_iq_settled and fault_v_settled
        if analysis_map[f"WT{i}_Signals_Settled"]:
            analysis_map[f"WT{i}_diq_pu"] = fault_iq_settled_value - pre_fault_iq_settled_value
            analysis_map[f"WT{i}_Fault_iq_pu"] = fault_iq_settled_value
            analysis_map[f"WT{i}_Fault_V_pu"] = fault_v_settled_value

            input_pars = (iq_pos_pu, fault_start_time, pre_fault_iq_settled_value, fault_iq_settled_time,
                      fault_iq_settled_value)
            wtg_iq_rise_time = signal_analysis.calc_aemo_step_rise_time(*input_pars)
            wtg_iq_settle_time = signal_analysis.calc_aemo_step_settling_time(*input_pars)

            wtg_max_iq_rise_time = max(wtg_max_iq_rise_time, wtg_iq_rise_time)
            wtg_max_iq_settle_time = max(wtg_max_iq_settle_time, wtg_iq_settle_time)


        all_wtg_settled = all_wtg_settled and analysis_map[f"WT{i}_Signals_Settled"]

    analysis_map[f"wtg_max_fault_iq_rise_time"] = wtg_max_iq_rise_time
    analysis_map[f"wtg_max_fault_iq_settle_time"] = wtg_max_iq_settle_time



    # Compute Bess Iq Values .
        
    bess_max_iq_rise_time = 0
    bess_max_iq_settle_time = 0
    all_bess_settled = True
    for i in [1, 2]:
        iq_pos_pu = data[f'BESS_Iq_pu:{i}']
        v_pos_pu = data[f'BESS_MV_V_pu:{i}']
        pre_fault_iq_settled, pre_fault_iq_settled_value, _ = \
            signal_analysis.calc_settled_value(iq_pos_pu, fault_start_time, min_settle_tol, ipoc_pu_tol)
        fault_iq_settled, fault_iq_settled_value, fault_iq_settled_time = \
            signal_analysis.calc_settled_value(iq_pos_pu, fault_end_time, min_settle_tol, ipoc_pu_tol)
        fault_v_settled, fault_v_settled_value, _ = \
            signal_analysis.calc_settled_value(v_pos_pu, fault_end_time, min_settle_tol, vpoc_pu_tol)

        analysis_map[f"BESS{i}_Signals_Settled"] = pre_fault_iq_settled and fault_iq_settled and fault_v_settled
        if analysis_map[f"BESS{i}_Signals_Settled"]:
            analysis_map[f"BESS{i}_diq_pu"] = fault_iq_settled_value - pre_fault_iq_settled_value
            analysis_map[f"BESS{i}_Fault_iq_pu"] = fault_iq_settled_value
            analysis_map[f"BESS{i}_Fault_V_pu"] = fault_v_settled_value

            input_pars = (iq_pos_pu, fault_start_time, pre_fault_iq_settled_value, fault_iq_settled_time,
                      fault_iq_settled_value)
            bess_iq_rise_time = signal_analysis.calc_aemo_step_rise_time(*input_pars)
            bess_iq_settle_time = signal_analysis.calc_aemo_step_settling_time(*input_pars)

            bess_max_iq_rise_time = max(bess_max_iq_rise_time, bess_iq_rise_time)
            bess_max_iq_settle_time = max(bess_max_iq_settle_time, bess_iq_settle_time)


        all_bess_settled = all_bess_settled and analysis_map[f"BESS{i}_Signals_Settled"]
        
    
    analysis_map[f"BESS_max_fault_iq_rise_time"] = bess_max_iq_rise_time
    analysis_map[f"BESS_max_fault_iq_settle_time"] = bess_max_iq_settle_time

           




    # Compute POC Iq Values.
    poc_iq_pos_pu = data['POC_Iq_pos_pu']
    # poc_i_pu = data['POC_I_kA'] / float(tuning['SYS_Ibase_kA'])

    # fault_ipoc_settled, fault_ipoc_settled_value, _ = \
    #     signal_analysis.calc_settled_value(poc_i_pu, fault_end_time, min_settle_tol, ipoc_pu_tol)
    # valid_test = valid_test and fault_ipoc_settled

    # if not fault_ipoc_settled:
    #     fault_ipoc_settled_value = poc_i_pu[fault_start_time:fault_end_time].values[-1]

    # analysis_map["Fault_Ipoc_mcc_pu"] = fault_ipoc_settled_value

    pre_fault_iq_settled, pre_fault_iq_settled_value, _ = \
        signal_analysis.calc_settled_value(poc_iq_pos_pu, fault_start_time, min_settle_tol, ipoc_pu_tol)
    fault_iq_settled, fault_iq_settled_value, fault_iq_settled_time = \
        signal_analysis.calc_settled_value(poc_iq_pos_pu, fault_end_time, min_settle_tol, ipoc_pu_tol)

    if not fault_iq_settled:
        fault_iq_settled_value = poc_iq_pos_pu[fault_start_time:fault_end_time].values[-1]


    valid_test = valid_test and pre_fault_iq_settled and fault_iq_settled and all_wtg_settled and all_bess_settled

    analysis_map["Prefault_iq_settled"] = pre_fault_iq_settled
    analysis_map["Prefault_iq_pu"] = pre_fault_iq_settled_value

    analysis_map["Fault_iq_settled"] = fault_iq_settled
    analysis_map["Fault_iq_pu"] = fault_iq_settled_value
    analysis_map["Fault_iq_mcc_pu"] = fault_iq_settled_value
    analysis_map["Delta_iq_mcc_pu"] = (fault_iq_settled_value - pre_fault_iq_settled_value)

    # Compute AEMO Rise and Settling time POC Iq if the signal was settled.
    if (pre_fault_iq_settled and fault_iq_settled) and (fault_start_time < fault_iq_settled_time):
        input_pars = (poc_iq_pos_pu, fault_start_time, pre_fault_iq_settled_value, fault_iq_settled_time,
                      fault_iq_settled_value)
        # print("calculating poc iq rise time")    
        iq_rise_time = signal_analysis.calc_aemo_step_rise_time(*input_pars)
        # print(f"iq rise is {iq_rise_time}")
        iq_settle_time = signal_analysis.calc_aemo_step_settling_time(*input_pars)

        mas_adeq_damped_diff, aas_adeq_damp_diff = signal_analysis.calc_aemo_adequately_damped(*input_pars)
        analysis_map["Fault_iq_AEMO_Rise_Time"] = iq_rise_time
        # Remove the 20ms 'rms filtering' of measurement delay from the settling time.
        # analysis_map["Fault_iq_AEMO_Settling_Time"] = iq_settle_time - (detected_fault_start_time - fault_start_time)
        analysis_map["Fault_iq_AEMO_Settling_Time"] = iq_settle_time
        analysis_map["Fault_iq_AEMO_MAS_Adeq_Damp"] = mas_adeq_damped_diff >= 0
        analysis_map["Fault_iq_AEMO_AAS_Adeq_Damp"] = aas_adeq_damp_diff >= 0


    vpoc_pu = data["POC_Vrms_pu"]
    
# #   Compute sum of all PCS RMS and check if it is less than 190% rating of battery current
#     PCS_I_BASE = 58.00297803
#     BESS_I_RATING = 48*1207
#     BESS_IRMS_MAX_TOTAL_PU = fault_data["BESS1_TERM_IRMS_PU"].max() + fault_data["BESS2_TERM_IRMS_PU"].max() + fault_data["BESS3_TERM_IRMS_PU"].max() + fault_data["BESS4_TERM_IRMS_PU"].max()
    
#     analysis_map["BESS_IRMS_MAX"] = BESS_IRMS_MAX_TOTAL_PU * PCS_I_BASE
#     analysis_map["BESS_Irms_over_190"] = analysis_map["BESS_IRMS_MAX"] > BESS_I_RATING*4*1.9
    

   # Compute POC Voltage Values.
    pre_fault_vpoc_settled, pre_fault_vpoc_settled_value, _ = \
        signal_analysis.calc_settled_value(vpoc_pu, fault_start_time, min_settle_tol, vpoc_pu_tol)
    fault_vpoc_settled, fault_vpoc_settled_value, _ = \
        signal_analysis.calc_settled_value(vpoc_pu, fault_end_time, min_settle_tol, vpoc_pu_tol)

    if not fault_vpoc_settled:
        fault_vpoc_settled_value = vpoc_pu[fault_start_time:fault_end_time].values[-1]

    valid_test = valid_test and pre_fault_vpoc_settled and fault_vpoc_settled
    analysis_map["Prefault_Vpoc_settled"] = pre_fault_vpoc_settled
    analysis_map["Prefault_Vpoc_pu"] = pre_fault_vpoc_settled_value
    analysis_map["Fault_Vpoc_settled"] = fault_vpoc_settled
    analysis_map["Fault_Vpoc_pu"] = fault_vpoc_settled_value
    analysis_map["Delta_V_pu(rel. 80%)"] = 0.80 - fault_vpoc_settled_value
    analysis_map["diq_mcc/dV"] = analysis_map["Delta_iq_mcc_pu"] / analysis_map["Delta_V_pu(rel. 80%)"]

    analysis_map["Fault_Vpoc_below_80"]= 0.8 > max(vpoc_pu[fault_start_time+0.03:fault_end_time])
    analysis_map["Fault_Vpoc_below_85"] = 0.85 > max(vpoc_pu[fault_start_time + 0.03:fault_end_time])

    # diff = vpoc_pu - data['WT_max_Vrms_All_Phases_pu']

    # analysis_map["POC_WTG_min_Diff"] = min(diff)
    # analysis_map["min_WT_max_Vrms_All_Phases_pu"] = min(data['WT_max_Vrms_All_Phases_pu'])


    # Compute Post Fault Active-Power Recovery Time.
    # Find the time it takes for the signal POC_P_MW to remain above 0.95*initial_ppoc_mw

    # Init. Ppoc.
    def calculate_95pc_recovery(sig, sig_recovery_val, sig_recovery_start_time, sig_recovery_end_time):
        if sig_recovery_val > 0:
            within_95pc_check = sig > 0.95 * sig_recovery_val
        else:
            within_95pc_check = sig < 0.95 * sig_recovery_val
        if list(sig[within_95pc_check].index) == []:
            first_recovery_time = None
        else:
            first_recovery_time = min(list(sig[within_95pc_check].index)) - sig_recovery_start_time

        within_95pc_check = within_95pc_check[::-1].iteritems()
        within_95pc_time = sig_recovery_end_time
        for curr_time, currently_above in within_95pc_check:
            if currently_above:
                within_95pc_time = curr_time
            else:
                break
        
        # if first_recovery_time is not None:
            # print(f"the active power recovery time is {first_recovery_time}")
        return first_recovery_time, (within_95pc_time - sig_recovery_start_time)
    # plt.plot(data["POC_P_MW"])
    initial_ppoc_mw = data["POC_P_MW"].values[0]
    sig = data["POC_P_MW"][post_fault_start_time:post_fault_end_time]
    first_recovery_time, full_recovery_time = calculate_95pc_recovery(sig, initial_ppoc_mw, post_fault_start_time,
                                                                      post_fault_end_time)
    # print(f"Ppoc_95pc_Recovery_Time is {full_recovery_time}")
    # plt.vlines(post_fault_start_time,0,1100,colors='r')
    # plt.vlines(post_fault_start_time+first_recovery_time,0,1100,colors='b')
    # plt.vlines(post_fault_start_time+full_recovery_time,0,1100,colors='k')
    # plt.hlines(initial_ppoc_mw*0.95,5.5,10,linestyles='dashed',colors='c')
    # plt.show()
    analysis_map["Ppoc_First_95pc_Recovery_Time"] = first_recovery_time
    analysis_map["Ppoc_95pc_Recovery_Time"] = full_recovery_time

    ppoc_first_within20MW_recovery_time, ppoc_within20MW_recovery_time = calculate_within_mw_recovery(sig, initial_ppoc_mw, post_fault_start_time,
                                                                      post_fault_end_time, 20)
    
    analysis_map[f"Ppoc_First_Within20MW_Recovery_Time"] = ppoc_first_within20MW_recovery_time
    analysis_map[f"Ppoc_Within20MW_Recovery_Time"] = ppoc_within20MW_recovery_time

    




    # Get Terminal Active-Power Recovery Times. 
    wtg_first_recovery_time = 0
    wtg_full_recovery_time = 0
    for sig_name in [f'WTG_MV_P_MW:{x}' for x in [1, 2, 3, 4, 5, 6]]:
        sig = data[sig_name][post_fault_start_time:post_fault_end_time]
        full_sig = data[sig_name]
        init_val = data[sig_name][fault_start_time:post_fault_end_time].values[0]
        first_recovery_time, full_recovery_time = calculate_95pc_recovery(sig, init_val, post_fault_start_time,
                                                                          post_fault_end_time)
        # print(f'wtg_first_recovery_time: {wtg_first_recovery_time}')
        # print(f'first_recovery_time: {first_recovery_time}')
        
        if not first_recovery_time or not full_recovery_time:
            # print("None!!")
            continue
        wtg_first_recovery_time = max(wtg_first_recovery_time, first_recovery_time)
        wtg_full_recovery_time = max(wtg_full_recovery_time, full_recovery_time)


        pre_fault_p_settled, pre_fault_p_settled_value, _ = \
            signal_analysis.calc_settled_value(full_sig, fault_start_time, min_settle_tol, 5) 
        
        fault_p_settled, fault_p_settled_value, fault_p_settled_time = \
            signal_analysis.calc_settled_value(full_sig, fault_end_time, min_settle_tol, 5)
        
        if fault_p_settled and pre_fault_p_settled:

            input_pars = (full_sig, fault_start_time, pre_fault_p_settled_value, fault_p_settled_time,
                fault_p_settled_value)

    analysis_map["WTGs_First_95pc_Recovery_Time"] = wtg_first_recovery_time
    analysis_map["WTGS_95pc_Recovery_Time"] = wtg_full_recovery_time


    # Get Terminal Active-Power Recovery Times 
    bess_first_recovery_time = 0
    bess_full_recovery_time = 0
    for sig_name in [f'BESS_MV_P_MW:{x}' for x in [1, 2]]:
        sig = data[sig_name][post_fault_start_time:post_fault_end_time]
        full_sig = data[sig_name]
        init_val = data[sig_name][fault_start_time:post_fault_end_time].values[0]
        first_recovery_time, full_recovery_time = calculate_95pc_recovery(sig, init_val, post_fault_start_time,
                                                                          post_fault_end_time)
        # print(f'bess_first_recovery_time: {bess_first_recovery_time}')
        # print(f'first_recovery_time: {first_recovery_time}')
        if not first_recovery_time or not full_recovery_time:
            # print("None!!")
            continue
        bess_first_recovery_time = max(bess_first_recovery_time, first_recovery_time)
        bess_full_recovery_time = max(bess_full_recovery_time, full_recovery_time)


        pre_fault_p_settled, pre_fault_p_settled_value, _ = \
            signal_analysis.calc_settled_value(full_sig, fault_start_time, min_settle_tol, 5) 
        
        fault_p_settled, fault_p_settled_value, fault_p_settled_time = \
            signal_analysis.calc_settled_value(full_sig, fault_end_time, min_settle_tol, 5)
        
        if fault_p_settled and pre_fault_p_settled:

            input_pars = (full_sig, fault_start_time, pre_fault_p_settled_value, fault_p_settled_time,
                fault_p_settled_value)

    analysis_map["BESS_First_95pc_Recovery_Time"] = bess_first_recovery_time
    analysis_map["BESS_95pc_Recovery_Time"] = bess_full_recovery_time




    # Iq AAS and MAS Settle-time. Not real settle time. Time till we are 'above' the mas/aas standard.
    mas_injection = max(min(2*(0.80 - fault_vpoc_settled_value), 1), 0)
    aas_injection = max(min(4 * (0.85 - fault_vpoc_settled_value), 1), 0)
    max_cont_curr = 1
    mas_iq_target = min(mas_injection + pre_fault_iq_settled_value, max_cont_curr)
    aas_iq_target = min(aas_injection + pre_fault_iq_settled_value, max_cont_curr)
    sig = poc_iq_pos_pu[fault_start_time:fault_end_time]

    above_mas_check = (sig > mas_iq_target)[::-1].iteritems()
    above_mas_time = fault_end_time
    for curr_time, currently_above in above_mas_check:
        if currently_above:
            above_mas_time = curr_time
        else:
            break
    analysis_map["Iq_Above_MAS_Time"] = above_mas_time - fault_start_time

    above_aas_check = (sig > aas_iq_target)[::-1].iteritems()
    above_aas_time = fault_end_time
    for curr_time, currently_above in above_aas_check:
        if currently_above:
            above_aas_time = curr_time
        else:
            break
    analysis_map["Iq_Above_AAS_Time"] = above_aas_time - fault_start_time

    # ------

    # # Maximum Compute Healthy Voltage on unaffected phases.
    # during_fault_healthy_map = healthy_voltage_map(data[fault_start_time:fault_end_time])
    # # post_fault_healthy_map = healthy_voltage_map(data[post_fault_start_time:post_fault_end_time])

    # analysis_map["Fault_Max_Healthy_Phase_pu"] = during_fault_healthy_map[fault_type]["max_healthy_pu"]
    # analysis_map["Fault_Min_Faulted_Phase_pu"] = during_fault_healthy_map[fault_type]["min_fault_pu"]
    # analysis_map["Post_Fault_Max_V_pu"] = max(data['POC_Vrms_pu'][post_fault_start_time:post_fault_end_time])

    # Check if the windfarm tripped.
    analysis_map["Valid_Test"] = valid_test

    # # Look over all FRT engages. and calculate lowest voltage.
    # for frt_sig in ['WT1_FRT', 'WT2_FRT', 'WT3_FRT', 'WT4_FRT']:
    #     lvrt_only = data[frt_sig].copy(deep=True)
    #     lvrt_only[data[frt_sig] == 1] = 0
    #     # lvrt_only = data[frt_sig][(data[frt_sig] == 0) | (data[frt_sig] == -1)]
    #     # hvrt_only = data[frt_sig][(data[frt_sig] == 0) | (data[frt_sig] == 1)]

    #     if frt_sig == 'WT1_FRT':
    #         lvrt_on_filter = (lvrt_only.diff() == -1)
    #         lvrt_off_filter = (lvrt_only.diff() == 1)
    #     else:
    #         lvrt_on_filter = lvrt_on_filter | (lvrt_only.diff() == -1)
    #         lvrt_off_filter = lvrt_off_filter | (lvrt_only.diff() == 1)

    # if len(data['POC_Vrms_pu'][lvrt_on_filter]) > 0:
    #     analysis_map["Min. WT_V_pu(LVRT_on)"] = min(data['WT1_V_pu'][lvrt_on_filter].values)
    #     analysis_map["Min. POC_Vrms_PU(LVRT_on)"] = min(data['POC_Vrms_pu'][lvrt_on_filter].values)
    # if len(data['POC_Vrms_pu'][lvrt_off_filter]) > 0:
    #     analysis_map["Min. WT_V_pu(LVRT_off)"] = min(data['WT1_V_pu'][lvrt_on_filter].values)
    #     analysis_map["Min. POC_Vrms_PU(LVRT_off)"] = min(data['POC_Vrms_pu'][lvrt_off_filter].values)

    # # The following code determines the accumulated time during the fault, that all 4 BESS governors were frozen.
    # fault_data = data[fault_start_time:fault_end_time]
    # time_delta = data.index[1] - data.index[0]
    # gov_frz_signals = [f'BESS{x}_GOV_FRZ' for x in [1, 2, 3, 4]]
    # fault_data['__Max_All_Gov_FRZ'] = fault_data[gov_frz_signals].max(axis=1)
    # fault_data['__Gov_FRZ_StintId'] = (fault_data['__Max_All_Gov_FRZ'].diff().fillna(0) != 0).cumsum()
    # # # cumsum() will create a 'stint ID' of sorts, so we can groupby it

    # stints = fault_data.groupby(['__Gov_FRZ_StintId', '__Max_All_Gov_FRZ'])['__Max_All_Gov_FRZ'].size()
    # stints = stints.loc[stints.index.get_level_values('__Max_All_Gov_FRZ') == 1].cumsum() * time_delta

    # analysis_map["Acc. GovFrz during fault"] = 0
    # if len(stints) > 0:
    #     analysis_map["Acc. GovFrz during fault"] = stints.values[0]

    # fault_data['__Max_VMP_FRZ'] = fault_data[[f'VMP_FRT']].max(axis=1)
    # fault_data['__VMP_FRZ_StintId'] = (fault_data['__Max_VMP_FRZ'].diff().fillna(0) != 0).cumsum()
    # # # cumsum() will create a 'stint ID' of sorts, so we can groupby it

    # stints = fault_data.groupby(['__VMP_FRZ_StintId', '__Max_VMP_FRZ'])['__Max_VMP_FRZ'].size()
    # stints = stints.loc[stints.index.get_level_values('__Max_VMP_FRZ') == 1].cumsum() * time_delta

    # analysis_map["Acc. VMPFrz during fault"] = 0
    # if len(stints) > 0:
    #     analysis_map["Acc. VMPFrz during fault"] = stints.values[0]


    # # Minimum Voltage during fault, and settled voltage during fault.
    # # This helps with discussing rise/settle-times of iq.
    # sig = data['POC_Vrms_pu'][detected_fault_start_time:fault_end_time]
    # analysis_map["Fault_Vpoc_Range_pu"] = max(sig) - min(sig)

    return analysis_map


def calculate_within_mw_recovery(sig, sig_recovery_val, sig_recovery_start_time, sig_recovery_end_time, threshold_MW):

        # PPOC Within MW threshold
        if sig_recovery_val > 0:
            within_mw_check = (sig > (sig_recovery_val - threshold_MW))
        else:
            within_mw_check = (sig < (sig_recovery_val + threshold_MW))

        if list(sig[within_mw_check].index) == []:
            sig_first_Within_MW_Recovery_Time = None
        else:
            # print(list(sig[within_mw_check].index))
            sig_first_Within_MW_Recovery_Time = \
                min(list(sig[within_mw_check].index)) - sig_recovery_start_time
        within_mw_check = within_mw_check[::-1].iteritems()
        within_mw_time = sig_recovery_end_time
        for curr_time, currently_above in within_mw_check:
            if currently_above:
                within_mw_time = curr_time
            else:
                break
        sig_within_MW_Recovery_Time = within_mw_time - sig_recovery_start_time

        return sig_first_Within_MW_Recovery_Time, sig_within_MW_Recovery_Time

if __name__ == '__main__':
    analysis_df = standard_analysis_main(
        file_name_prefix=["s5255_Fault_Ures", "s5255_Fault_Ohms"], per_scenario_analysis=per_scenario_analysis,
        summary_analysis=None)

