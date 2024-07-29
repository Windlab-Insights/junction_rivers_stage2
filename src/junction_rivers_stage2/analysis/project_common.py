import numpy as np
from  junction_rivers.analysis import signal_analysis


def check_all_tx_settled(data):
    tx_settled = []
    for idx in range(1, 7):
        tx_settled.append(signal_analysis.signal_is_constant(data[f'TX{idx}_Settled'], value=1))
    return all(tx_settled)


def check_wtg_no_trip(data):
    wtg_no_trip = []
    for idx in range(1, 7):
        wtg_no_trip.append(signal_analysis.signal_is_constant(data[f'WTG_Trip_Flag:1:{idx}'], value=0))
    return all(wtg_no_trip)


def check_bess_no_trip(data):
    bess_no_trip = []
    for idx in range(1, 3):
        bess_no_trip.append(signal_analysis.signal_is_constant(data[f'BESS_Trip_Flag:1:{idx}'], value=0))
    return all(bess_no_trip)


def check_wtg_no_hvrt(data):
    wtg_no_frt = []
    for idx in range(1, 7):
        wtg_no_frt.append(signal_analysis.signal_is_constant(data[f'WTG_HVRT_Flag:1:{idx}'], value=0))
    return all(wtg_no_frt)

def check_wtg_no_lvrt(data):
    wtg_no_frt = []
    for idx in range(1, 7):
        wtg_no_frt.append(signal_analysis.signal_is_constant(data[f'WTG_LVRT_Flag:1:{idx}'], value=0))
    return all(wtg_no_frt)


def check_bess_no_frt(data):
    bess_no_frt = []
    for idx in range(1, 3):
        bess_no_frt.append(signal_analysis.signal_is_constant(data[f'BESS{idx}_AVR_FRZ'], value=0))
        bess_no_frt.append(signal_analysis.signal_is_constant(data[f'BESS{idx}_GOV_FRZ'], value=0))
    return all(bess_no_frt)


def check_vmp_no_frt(data):
    return signal_analysis.signal_is_constant(data[f'VMP_FRT'], value=0)

def qref_mvar_via_voltage_droop(tuning, vref_pu, vpoc_pu):
    required_tuning_keys = ['SYS_Qbase_MVAr', 'SYS_Qmin_MVAr', 'SYS_Qmax_MVAr', 'SYS_Vdroop_perc']

    if not all([x in tuning for x in required_tuning_keys]):
        raise ValueError(f'Tuning File must contain: {required_tuning_keys}')

    qbase_mvar = float(tuning['SYS_Qbase_MVAr'])
    qmin_mvar = float(tuning['SYS_Qmin_MVAr'])
    qmax_mvar = float(tuning['SYS_Qmax_MVAr'])
    vdroop_ratio = float(tuning['SYS_Vdroop_perc'])
    return np.clip((1/vdroop_ratio)*(vref_pu - vpoc_pu)*qbase_mvar, qmin_mvar, qmax_mvar)


def adjusted_pref_given_grid_hz(tuning, pref_wtg_mw, pref_bess_mw, grid_hz):
    required_tuning_keys = ['BESS_PPC_Pmax_MW', 'BESS_PPC_Pmin_MW', 'SYS_Pbase_MW', 'VSG_F_Deadband_Hz', 'SYS_Fbase_Hz',
                            'SYS_F_Droop_perc']

    if not all([x in tuning for x in required_tuning_keys]):
        raise ValueError(f'Tuning File must contain: {required_tuning_keys}')

    CAP_BESS_MAX_MW = float(tuning['BESS_PPC_Pmax_MW']) * 4
    CAP_BESS_MIN_MW = float(tuning['BESS_PPC_Pmin_MW']) * 4
    CAP_EMP_MAX_MW = float(tuning['SYS_Pbase_MW'])
    CAP_POC_MAX_MW = CAP_EMP_MAX_MW
    CAP_POC_MIN_MW = CAP_BESS_MIN_MW
    FBASE_HZ = float(tuning['SYS_Fbase_Hz'])
    VSG_F_DEADBAND_PU = float(tuning['VSG_F_Deadband_Hz']) / FBASE_HZ
    BESS_F_DEADBAND_PU = VSG_F_DEADBAND_PU

    # Notice, we are splitting the droop over two units, so we multiply by 2.
    fdroop_on_pbase = float(tuning['SYS_F_Droop_perc']) * 2

    delta_f_pu = (grid_hz - FBASE_HZ) / FBASE_HZ

    pref_total = pref_wtg_mw + pref_bess_mw

    pos_rem_mw = CAP_POC_MAX_MW - pref_total
    neg_rem_mw = CAP_POC_MIN_MW - pref_total

    total_bess_max_mw = min(pref_bess_mw + pos_rem_mw, CAP_BESS_MAX_MW)
    total_bess_min_mw = max(pref_bess_mw + neg_rem_mw, CAP_BESS_MIN_MW)

    total_wtg_max_mw = CAP_POC_MAX_MW - total_bess_max_mw
    total_wtg_min_mw = 0

    vsg_bess_max_mw = total_bess_max_mw - pref_bess_mw
    vsg_bess_min_mw = total_bess_min_mw - pref_bess_mw

    vsg_wtg_max_mw = total_wtg_max_mw - pref_wtg_mw
    vsg_wtg_min_mw = total_wtg_min_mw - pref_wtg_mw

    wtg_droop_norm_on_pbase = fdroop_on_pbase * (CAP_EMP_MAX_MW / CAP_POC_MAX_MW)
    bess_droop_norm_on_pbase = fdroop_on_pbase * (CAP_BESS_MAX_MW / CAP_POC_MAX_MW)

    pos_adjusted_wtg_droop = wtg_droop_norm_on_pbase
    neg_adjusted_wtg_droop = wtg_droop_norm_on_pbase
    pos_adjusted_bess_droop = bess_droop_norm_on_pbase
    neg_adjusted_bess_droop = bess_droop_norm_on_pbase

    def calc_fdroop_delta_p(delta_f, pos_droop_ratio, neg_droop_ratio, dead_band_pu, pmax, pref, vsg_min, vsg_max):
        if abs(delta_f) <= dead_band_pu:
            # Within Deadband.
            delta_p = 0
        else:
            # Not within deadband.
            if delta_f >= 0:
                # positive delta_f. Compute value, we will saturate at the end.
                delta_p = -(pmax / neg_droop_ratio) * delta_f
            else:
                # negative delta_f
                delta_p = -(pmax / pos_droop_ratio) * delta_f

        delta_p = np.clip(delta_p, vsg_min, vsg_max)
        return delta_p

    wtg_droop_delta_p = calc_fdroop_delta_p(delta_f_pu, pos_adjusted_wtg_droop, neg_adjusted_wtg_droop, VSG_F_DEADBAND_PU,
                                            CAP_EMP_MAX_MW, pref_wtg_mw, vsg_wtg_min_mw, vsg_wtg_max_mw)
    bess_droop_delta_p = calc_fdroop_delta_p(delta_f_pu, pos_adjusted_bess_droop, neg_adjusted_bess_droop,
                                             BESS_F_DEADBAND_PU, CAP_BESS_MAX_MW, pref_bess_mw, vsg_bess_min_mw,
                                             vsg_bess_max_mw)
    poc_delta_p = wtg_droop_delta_p + bess_droop_delta_p
    adj_poc_pref = poc_delta_p + pref_bess_mw + pref_wtg_mw
    adj_wtg_pref = wtg_droop_delta_p + pref_wtg_mw
    adj_bess_pref = bess_droop_delta_p + pref_bess_mw
    return adj_poc_pref, adj_wtg_pref, adj_bess_pref
