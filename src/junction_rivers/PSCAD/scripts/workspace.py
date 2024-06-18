from collections import OrderedDict
import pandas as pd
from junction_rivers.Analysis import Analysis
import os
import numpy as np
from typing import List, Optional, Tuple, Union, Dict
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as backend_pdf
import matplotlib.gridspec as gridspec

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
    
    pdf_path = "C:\Temp\plot_test.pdf"
    pdf_path1 = "C:\Temp\plot_test1.pdf"
    
    fig1 = plt.figure()
    fig2 = plt.figure()

    x = np.array([[1,4],[1,2],[1,9]])
    y = np.array([[1,2],[1,2],[1,2]])
    
    layout_main_1 = gridspec.GridSpec(1, 1, figure=fig1)
    layout_main_2 = gridspec.GridSpec(1, 1, figure=fig2)
    column_1 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=layout_main_1[0])
    column_2 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=layout_main_2[0])
    
    ax_1: List[plt.Axes] = []
    ax_2: List[plt.Axes] = []

    ax_1.append(fig1.add_subplot(column_1[0]))
    ax_2.append(fig2.add_subplot(column_2[0]))
    
    ax_1[0].plot(y)
    ax_2[0].plot(x)
    
    fig1.savefig(pdf_path1)
    pdf = backend_pdf.PdfPages(pdf_path)
    try:
        with backend_pdf.PdfPages(pdf_path) as pdf:
            pdf.savefig(fig1)
            pdf.savefig(fig2)
            print("### Plot saved to: "+ pdf_path)
    except Exception as e:
        print(f"### pdf save didn't work {e}")
        
    plt.cla() 
    fig1.clf()
    fig2.clf()
    plt.close()
    
    
