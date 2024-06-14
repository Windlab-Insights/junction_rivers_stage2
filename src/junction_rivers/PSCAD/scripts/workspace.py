from collections import OrderedDict
import pandas as pd
from junction_rivers.Analysis import Analysis
import os

def workspace():
    df = pd.DataFrame()
    Dict = OrderedDict()
    Dict['Title A'] = "aaa"
    Dict['Title B'] = "bbb"
    print(f" orderedDict: {Dict}")
    df = df.append(other=pd.DataFrame([Dict]), ignore_index=True)
    print(f"df: {df}")
    Dict2 = OrderedDict()
    Dict2['Title A'] = "aab"
    Dict2['Title B'] = "bbc"
    print(f" orderedDict2: {Dict2}")
    df = df.append(other=pd.DataFrame([Dict2]), ignore_index=True)
    print(f"df: {df}")
    print(f"df[A]: {df['Title A']}")
    data (df)
    print(f"df: {df}")
    
def data(df: pd.DataFrame):
    Dict3 = OrderedDict()
    Dict3['Title A'] = "acc"
    Dict3['Title B'] = "bcd"
    print(f" orderedDict3: {Dict3}")
    df = df.append(other=pd.DataFrame([Dict3]), ignore_index=True)
    print(f"df: {df}")
    
if __name__ == '__main__':
    path = os.path.join("C:\Temp","fdroop_plotter_test.csv")
    analysis = Analysis()
    analysis_map = OrderedDict()
    analysis_map['Freq_During_Dist_Hz'] = 1
    analysis_map['POC_P_During_Dist_MW'] = 2
    analysis_map['WTG_P_During_Dist_MW'] = 3
    analysis_map['BESS_P_During_Dist_MW'] = 4

    analysis_map['Delta_Freq_During_Dist_Hz'] = 5
    analysis_map['Delta_POC_P_During_Dist_MW'] = 6
    analysis_map['Delta_WTG_P_During_Dist_MW'] = 7
    analysis_map['Delta_BESS_P_During_Dist_MW'] = 8
    df = pd.DataFrame([analysis_map])
    analysis_map2 = OrderedDict()
    analysis_map2['Freq_During_Dist_Hz'] = 12
    analysis_map2['POC_P_During_Dist_MW'] = 22
    analysis_map2['WTG_P_During_Dist_MW'] = 32
    analysis_map2['BESS_P_During_Dist_MW'] = 42

    analysis_map2['Delta_Freq_During_Dist_Hz'] = 52
    analysis_map2['Delta_POC_P_During_Dist_MW'] = 62
    analysis_map2['Delta_WTG_P_During_Dist_MW'] = 72
    analysis_map2['Delta_BESS_P_During_Dist_MW'] = 82
    df = df.append(pd.DataFrame([analysis_map2]))
    analysis.set_fdroop_log(df)
    analysis.plot_fdroop(path)
    
