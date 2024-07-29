import math
import numpy as np
import pandas as pd
from scipy.signal import butter, lfilter
import matplotlib.pyplot as plt


def signal_is_constant(raw_signal: pd.Series, value=None, variation_tolerance=0, value_tolerance=0):
    is_constant = False
    total_variation = abs(raw_signal.max() - raw_signal.min())
    if total_variation <= variation_tolerance:
        is_constant = True

    # Reversing logic so that signal_is_constant can be used to detect constancy only, without needing to
    # provide the value.
    has_value = True
    if value is not None and abs(value - raw_signal.mean()) > value_tolerance:
        has_value = False

    return is_constant and has_value


# def calc_natural_frequency(signal: pd.Series, sampling_rate):
#     # Compute power spectral density
#     psd_freq, psd_pxx = scipy.signal.welch(signal, fs=sampling_rate)
#     pxx = pd.Series(psd_pxx, index=psd_freq)
#     natural_frequency = pxx.argmax()
#     return natural_frequency

def calc_settled_value(raw_signal: pd.Series, settledby_time: float, time_tol: float, val_tol: float):
    """
    Args:
        signal: The input pd.Series signal.
        settledby_time: The time before which the settling value is to be calculated.
        time_tol: The minimum time that a signals variation must be less than val_tol before being concluded as settled.
        val_tol: The maximum acceptable variation for a signal to be concluded as settled.

    Returns:
        signal_settled: Returns True of signal is settled, otherwise False.
         settled_value: Returns the value of the settled signal if it settled.
         settled_time: Returns the earliest time that the signal remains within val_tol of settled_valued, if settled.
    """

    signal = raw_signal[:settledby_time].iloc[:-1]

    settled_value = -1
    settling_time = -1.
    signal_settled = False
    if time_tol < 0 or val_tol < 0:
        raise ValueError(f"Tolerances must be positive values. Supplied (xtol, ytol)=({time_tol},{val_tol})")
    else:
        # In the case that val_tol is 0, increase it so that the floating point comparisons work.
        val_tol = 0.0000001 if val_tol == 0 else val_tol
        settle_start = settledby_time - time_tol
        # Very last point has been removed as it may be contminated with the following step signal
        xtruncated = signal[settle_start:settledby_time]
        # print(settle_start,settledby_time)
        diff = max(xtruncated) - min(xtruncated)
        if diff < val_tol or math.isclose(diff, val_tol):
            signal_settled = True
            settled_value = signal[:settledby_time].iloc[-1]
            upper_settling_boundary = settled_value + (val_tol / 2)
            lower_settling_boundary = settled_value - (val_tol / 2)
            settling_check = signal[:settledby_time].between(lower_settling_boundary,
                                                             upper_settling_boundary)[::-1].iteritems()
            for curr_time, currently_settled in settling_check:
                if currently_settled:
                    settling_time = curr_time
                else:
                    break

    return signal_settled, settled_value, settling_time

def low_pass_filter(signal: pd.Series, cut_off: float):
        fs = 1 / (signal.index[3] - signal.index[2])

        b, a = butter(1, cut_off, fs=fs, btype='low', analog=False)
        filtered_signal = lfilter(b, a, signal.values - signal.values[0]) + signal.values[0]

        return pd.Series(
            data=filtered_signal,
            index=signal.index,
        )

def calc_aemo_step_rise_time(signal: pd.Series, pre_step_settled_time: float, pre_step_settled_value: float,
                             post_step_settled_time: float, post_step_settled_value: float):
    # rise time
    # In relation to a control system, the time taken for an output quantity to rise from
    # 10% to 90% of the maximum change induced in that quantity by a step change of
    # an input quantity. This assumes that systems are non-minimal phase, i.e., they don't undershoot
    # before rising towards their settled value.

    # Calculate Maximum Change.
    rise_time =-1
    signal_window = signal[pre_step_settled_time:post_step_settled_time]

    if len(signal_window) == 0:
        return 0

    # max_peak_time = signal_window.idxmax()
    # max_peak_val = signal_window[max_peak_time]
    # min_peak_time = signal_window.idxmin()
    # min_peak_val = signal_window[min_peak_time]
    signal_window_filtered= low_pass_filter(signal_window, 50)
    max_over_window = max(signal_window_filtered)
    min_over_window = min(signal_window_filtered)
    maximum_change_induced = max_over_window - min_over_window

    if post_step_settled_value > pre_step_settled_value:
        # Stepping UP.
        milestone1 = min_over_window + (0.1 * maximum_change_induced)
        milestone2 = min_over_window + (0.9 * maximum_change_induced)

        last_time_below_milestone1 = signal_window[signal_window <= milestone1].index[0]
        first_time_above_milestone2 = signal_window[signal_window >= milestone2].index[0]

        # fig = plt.figure()  # subplot(number_of_rows, 4)
        # print(first_time_above_milestone2 - last_time_below_milestone1)
        # plt.plot(signal_window)
        # plt.plot(signal_window_filtered)
        # plt.vlines(last_time_below_milestone1,min_over_window,max_over_window,colors='r')
        # plt.vlines(first_time_above_milestone2,min_over_window,max_over_window,colors='b')
        # plt.show()
        rise_time = first_time_above_milestone2 - last_time_below_milestone1
        # print(f"in signal analysis rise time is {rise_time}")
        return rise_time
    else:
        # Stepping Down.
        milestone1 = max_over_window - (0.1 * maximum_change_induced)
        milestone2 = max_over_window - (0.9 * maximum_change_induced)
        last_time_above_milestone1 = signal_window[signal_window >= milestone1].index[0]
        first_time_below_milestone2 = signal_window[signal_window <= milestone2].index[0]
        return first_time_below_milestone2 - last_time_above_milestone1
    return -1


def calc_aemo_step_settling_time(signal: pd.Series, pre_step_settled_time: float, pre_step_settled_value: float,
                                 post_step_settled_time: float, post_step_settled_value: float):
    # AEMO Standard.
    # In relation to a control system, the time measured from initiation of a step change
    # in an input quantity to the time when the magnitude of error between the output
    # quantity and its final settling value remains less than 10% of:
    # (a) if the sustained change in the quantity is less than half of the maximum
    # change in that output quantity, the maximum change induced in that output
    # quantity; or
    # (b) the sustained change induced in that output quantity.

    # Calculate Maximum Change.
    signal_window = signal[pre_step_settled_time:post_step_settled_time]

    if len(signal_window) == 0:
        return 0

    max_over_window = max(signal_window)
    min_over_window = min(signal_window)
    maximum_change_induced = max_over_window - min_over_window
    sustained_change = abs(post_step_settled_value - pre_step_settled_value)
    if sustained_change < (maximum_change_induced / 2):
        aemo_change = maximum_change_induced
    else:
        aemo_change = sustained_change

    upper_settling_boundary = post_step_settled_value + (aemo_change * 0.1)
    lower_settling_boundary = post_step_settled_value - (aemo_change * 0.1)

    # Note: [::-1] has the effect of reversing the list.
    # What we are doing is checking which values in signal_window are between upper and lower boundary.
    # this then creates a Series of booleans. We then reverse the time-order using [::-1], and then
    # iterate backwards to find the first False. That False is the last time we are outside of settling
    # range.
    settled_check = signal_window.between(lower_settling_boundary, upper_settling_boundary)[::-1].iteritems()
    settled_time = post_step_settled_time
    for curr_time, currently_settled in settled_check:
        if currently_settled:
            settled_time = curr_time
        else:
            break
    # Finally, Settling time is the difference
    aemo_settling_time = settled_time - pre_step_settled_time
    return aemo_settling_time


def calc_transient_characteristics(signal: pd.Series, pre_step_settled_time: float,
                                   pre_step_settled_value: float,
                                   post_step_settled_time: float, post_step_settled_value: float):

    signal_window = signal[pre_step_settled_time:post_step_settled_time]

    if post_step_settled_value - pre_step_settled_value > 0:
        step_direction = 1
        time_of_peak = signal_window.idxmax()
        time_of_start_climb = signal_window.idxmin()

    else:
        step_direction = -1
        time_of_peak = signal_window.idxmin()
        time_of_start_climb = signal_window.idxmax()

    peak_value = signal_window[time_of_peak]

    if (peak_value < post_step_settled_value and step_direction == -1) or (peak_value > post_step_settled_value and step_direction == 1):
        # Underdamped.
        percentage_overshoot = abs((peak_value - post_step_settled_value) / post_step_settled_value)
        damped_frequency = np.pi / (time_of_peak - time_of_start_climb)
        damping_ratio = -np.log(percentage_overshoot) / np.sqrt(np.pi ** 2 + np.log(percentage_overshoot) ** 2)
        decay_rate = (damping_ratio / (np.sqrt(1 - damping_ratio ** 2))) * damped_frequency
        halving_time = np.log(2) / decay_rate
    else:
        # Overdamped
        decay_rate, damped_frequency, damping_ratio, percentage_overshoot, halving_time = (0, 0, 0, 0, 0)

    return decay_rate, damped_frequency, damping_ratio, percentage_overshoot, halving_time


def calc_aemo_adequately_damped_decay_rate(omega_rads):
    # Automatic Access Standard.
    # if omega_n in [0, 0.05)     min_aas_decay_rate = 0.43643578 * omega_n
    # if omega_n in [0.05, 0.6)   min_aas_decay_rate = 0.14
    # if omega_n in [0.6, infty)  min_aas_decay_rate = 0.100503782 * omega_n
    def min_aas_decay_rate(omega):
        result = 0
        if 0 <= omega < 2 * np.pi * 0.05:
            result = 0.43643578 * omega
        if 2 * np.pi * 0.05 <= omega < 2 * np.pi * 0.6:
            result = 0.14
        if 2 * np.pi * 0.6 <= omega:
            result = 0.100503782 * omega
        return result

    # Minimum Access Standard.
    # if omega_n in [0, 0.05)     min_mas_decay_rate = 0.43643578 * omega_n
    # if omega_n in [0.05, 0.6)   min_mas_decay_rate = 0.14
    # if omega_n in [0.6, infty)  min_mas_decay_rate = 0.050062617 * omega_n
    def min_mas_decay_rate(omega):
        result = 0
        if 0 <= omega < 2 * np.pi * 0.05:
            result = 0.43643578 * omega
        if 2 * np.pi * 0.05 <= omega < 2 * np.pi * 0.6:
            result = 0.14
        if 2 * np.pi * 0.6 <= omega:
            result = 0.050062617 * omega
        return result

    return min_mas_decay_rate(omega_rads), min_aas_decay_rate(omega_rads)


def calc_aemo_adequately_damped(signal: pd.Series, pre_step_settled_time: float, pre_step_settled_value: float,
                                post_step_settled_time: float, post_step_settled_value: float):
    """
    Given a step response, calculate if it passes either the minimum or automatic acceptance standard for being
    adequately damped.
    Args:
        signal: pd.Signal
        pre_step_settled_time: int
        pre_step_settled_value: int
        post_step_settled_time: int
        post_step_settled_value: int

    Returns:
        passed_minimum_acceptance_standard : bool
        passed_automatic_acceptance_standard : bool
    """
    decay_rate, damped_frequency, _, _, _ = calc_transient_characteristics(signal, pre_step_settled_time,
                                                                           pre_step_settled_value,
                                                                           post_step_settled_time,
                                                                           post_step_settled_value)
    min_mas_decay_rate, min_aas_decay_rate = calc_aemo_adequately_damped_decay_rate(damped_frequency)
    return decay_rate - min_mas_decay_rate, decay_rate - min_aas_decay_rate


def calc_aemo_step_characteristics(signal: pd.Series, pre_step_settled_time: float,
                                   post_step_settled_time: float, settle_time_tol: float, settle_val_tol: float):
    pre_settled, pre_settled_value, _ = calc_settled_value(signal, pre_step_settled_time, settle_time_tol, settle_val_tol)
    post_settled, post_settled_value, _ = calc_settled_value(signal, post_step_settled_time, settle_time_tol, settle_val_tol)

    if pre_settled and post_settled:
        rise_time = calc_aemo_step_rise_time(signal, pre_step_settled_time, pre_settled_value, post_step_settled_time,
                                             post_settled_value)
        settle_time = calc_aemo_step_settling_time(signal, pre_step_settled_time, pre_settled_value,
                                                   post_step_settled_time, post_settled_value)

        min_damping_standard, auto_damping_standard = calc_aemo_adequately_damped(signal,
                                                                                  pre_step_settled_time,
                                                                                  pre_settled_value,
                                                                                  post_step_settled_time,
                                                                                  post_settled_value)

        return pre_settled_value, post_settled_value, rise_time, settle_time, min_damping_standard, \
               auto_damping_standard
    else:
        return pre_settled_value, post_settled_value, -1, -1, -1, -1
    
    
# Expected Active Power Freq Droop  
def get_expected_fdroop_signal(spec_dict: dict, df: pd.DataFrame) -> pd.Series : 
    poc_pbase_mw = float(spec_dict["substitutions"]["SYS_Pbase_MW"])
    VSG_Fdroop_on_Pbase = float(spec_dict["substitutions"]["VSG_Fdroop_on_Pbase"])
    VSG_F_Deadband_Hz = float(spec_dict["substitutions"]["VSG_F_Deadband_Hz"])
    poc_pref_mw = df["POC_WTG_Pref_MW"]
    poc_freq_hz = df['POC_Freq_Hz']
    poc_ferr_hz = (poc_freq_hz - 50)
    poc_ferr_hz[poc_ferr_hz.abs()<VSG_F_Deadband_Hz] = 0
    signal = poc_pref_mw - poc_ferr_hz * VSG_Fdroop_on_Pbase * poc_pbase_mw
    return signal


# Expected Reactive Power Voltage Droop
def get_expected_vdroop_signal(spec_dict: dict, df: pd.DataFrame) -> pd.Series : 
    poc_vref_pu = df["POC_Vref_pu"]
    poc_vrms_pu = df["POC_Vrms_pu"]
    poc_qbase_mvar = float(spec_dict["substitutions"]["SYS_Qbase_MVAr"])
    SYS_Vdroop_on_Qbase = float(spec_dict["substitutions"]["SYS_Vdroop_on_Qbase"])
    poc_verr_pu = poc_vref_pu - poc_vrms_pu
    raw_Qref_droop_Ideal_MVAr = poc_vref_pu + (poc_verr_pu / SYS_Vdroop_on_Qbase) * poc_qbase_mvar
    raw_Qref_droop_Ideal_MVAr[raw_Qref_droop_Ideal_MVAr > poc_qbase_mvar] = poc_qbase_mvar
    raw_Qref_droop_Ideal_MVAr[raw_Qref_droop_Ideal_MVAr < -poc_qbase_mvar] = -poc_qbase_mvar
    signal = raw_Qref_droop_Ideal_MVAr
    return signal

    
    

# asdf
#     # # Empirically Confirm decay rate.
#     # abs_normalised_transient = abs(signal[time_to_peak:post_step_settled_time] - post_step_settled_value) / \
#     #                            abs(peak_value - post_step_settled_value)
#     #
#     # # Now tweak the decay_rate
#     # right_shift = 0.01 * (post_step_settled_time - time_to_peak)
#     # abs_settle_tol = 0.01 * post_step_settled_value
#     # emp_decay_rate = theoretical_decay_rate
#     # first_adjust = True
#     # for curr_time, curr_value in abs_normalised_transient.iteritems():
#     #     envelope = np.exp(-emp_decay_rate * (curr_time - time_to_peak - right_shift)) + abs_settle_tol
#     #     if curr_value - envelope > 0.01:
#     #         if first_adjust:
#     #             # Need to adjust envelope
#     #             abs_settle_tol = abs_settle_tol + curr_value - envelope
#     #             first_adjust = False
#     #         else:
#     #             emp_decay_rate = -np.log(curr_value - abs_settle_tol) / (curr_time - time_to_peak)
#     #
#     # # ---
#     # fig, axes = plt.subplots()
#     # abs_normalised_transient.plot(ax=axes)
#     # envelope = pd.Series(
#     #     np.exp(-emp_decay_rate * (abs_normalised_transient.index - time_to_peak)) + abs_settle_tol,
#     #     index=abs_normalised_transient.index)
#     # envelope.plot(ax=axes)
#     # plt.savefig("test.png")
#     # plt.show()
#     # # ---
#     #
#     # return emp_decay_rate, theoretical_decay_rate, damping_ratio, damped_frequency, percentage_overshoot
#     #
#     #
#     # # Now compute percentage overshoot, which gives theoretical damping ratio.
#     #
#     # # NOTE: At this moment, my understanding of halving time requires that a signal is in steady-state, it is then
#     # #       perturbed via some disturbance, and then recovers back to its steady state. The halving time is the
#     # #       interval starting at the moment deviation is the largest away from steady-state, and finishes at the moment
#     # #       that  the signal's deviation remains within 50% of that worst-case deviation.
#     #
#     # # Calculate the maximum deviation time and value.
#     # change_sign = np.sign(post_step_settled_value - pre_step_settled_value)
#     # signal_window = signal[pre_step_settled_time:post_step_settled_time]
#     # if max(signal_window)) > abs(post_step_settled_value - pre_step_settled_value):
#     #
#     #
#     # time_of_max_deviation = change_sign*signal_window.argmax()
#     #
#     #
#     #
#     # error_signal = abs(signal_window - post_step_settled_value)
#     # max_deviation_time = error_signal.argmax()
#     # max_deviation = error_signal[max_deviation_time]
#     #
#     # # Calculate the halving time.
#     # half_max_deviation = max_deviation / 2
#     #
#     # lower_boundary = steady_state_value - half_max_deviation
#     # upper_boundary = steady_state_value + half_max_deviation
#     # # Note: [::-1] makes the iterator loops in reverse order.
#     # settled_check = signal_window.between(lower_boundary, upper_boundary)[::-1].iteritems()
#     # settled_time = None
#     # for curr_time, currently_settled in settled_check:
#     #     if currently_settled:
#     #         settled_time = curr_time
#     #     else:
#     #         break
#     # if settled_time is None:
#     #     return -1
#     # else:
#     #     halving_time = settled_time - max_deviation_time
#     #     return halving_time
