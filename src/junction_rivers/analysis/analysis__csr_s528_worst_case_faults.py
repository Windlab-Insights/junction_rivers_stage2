import pandas as pd
from collections import OrderedDict
import signal_analysis
from standard_mains import standard_analysis_main, find_filepath_beneath_dir
from project_common import check_all_tx_settled, check_wtg_no_trip, check_bess_no_trip, check_wtg_no_frt, check_bess_no_frt, check_vmp_no_frt


"""
Low-Voltage Faults.
"""
#

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

        max_healthy_phase = max(max(data[f"POC_PhaseA_V_pu"]),
                                max(data[f"POC_PhaseB_V_pu"]),
                                max(data[f"POC_PhaseC_V_pu"]))
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

    model_init_time = float(tuning['TIME_Full_Init_Time_sec'])
    post_init_time = model_init_time + float(scenario['Post_Init_Duration_s'])

    fault_start_time = model_init_time + float(scenario['Fault_Time'])
    fault_end_time = fault_start_time + float(scenario['Fault_Duration'])
    post_fault_end_time = fault_end_time + float(scenario['Post_Fault_Duration'])

    data = run_data[model_init_time:post_init_time]
    analysis_map = OrderedDict()

    analysis_map['Whole_Run_Max_POC_PhaseA_Iinst_kA'] = max(data[f'POC_Iinst_kA:1'])
    analysis_map['Whole_Run_Max_POC_PhaseB_Iinst_kA'] = max(data[f'POC_Iinst_kA:2'])
    analysis_map['Whole_Run_Max_POC_PhaseC_Iinst_kA'] = max(data[f'POC_Iinst_kA:3'])
    analysis_map['Whole_Run_Max_POC_Iinst_kA'] = max(max(data[f'POC_Iinst_kA:1']),
                                           max(data[f'POC_Iinst_kA:2']),
                                           max(data[f'POC_Iinst_kA:3']))

    fault_data = data[fault_start_time:fault_end_time]
    analysis_map['Fault_Max_POC_PhaseA_Iinst_kA'] = max(fault_data[f'POC_Iinst_kA:1'])
    analysis_map['Fault_Run_Max_POC_PhaseB_Iinst_kA'] = max(fault_data[f'POC_Iinst_kA:2'])
    analysis_map['Fault_Run_Max_POC_PhaseC_Iinst_kA'] = max(fault_data[f'POC_Iinst_kA:3'])
    analysis_map['Fault_Run_Max_POC_Iinst_kA'] = max(max(fault_data[f'POC_Iinst_kA:1']),
                                           max(fault_data[f'POC_Iinst_kA:2']),
                                           max(fault_data[f'POC_Iinst_kA:3']))
    return analysis_map

# def summary_analysis(tuning, analysis_df, data_root_dir, output_dir_path):


if __name__ == '__main__':
    analysis_df = standard_analysis_main(
        file_name_prefix=["s528_Faults_Ures", "s528_Fault_Ohms", "s5255_Faults_Ures", "s5255_Fault_Ohms"], per_scenario_analysis=per_scenario_analysis,
        summary_analysis=None)

