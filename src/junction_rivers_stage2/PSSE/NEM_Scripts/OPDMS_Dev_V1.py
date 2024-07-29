#-------------------------------------------------------------------------------
# Name:        Base Case developement
# Purpose:
#
# Author:      Sajjad Hadavi
#
# Created:     210/10/2023
# Copyright:   (c) Sajjad Hadavi 2023
# Licence:     <your licence>


#-------------------------------------------------------------------------------

#%% Import Libraries:
# import pssepath
# pssepath.add_pssepath()
import argparse
import time
import datetime
import pandas as pd
import os, sys, pdb
import xlrd
import csv
import re
from string import digits
import matplotlib
import networkx as nx
from matplotlib.ticker import FormatStrFormatter
from matplotlib import ticker as mtick
from matplotlib.font_manager import FontProperties
from matplotlib.backends.backend_pdf import PdfPages as pdf
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np
import math
import psse34
import dyntools as dt
import matplotlib.image as image
import redirect
import psspy
from psspy import _i,_f,_s
import glob
import scipy.signal
import warnings
import numpy.polynomial.polynomial as poly
import shutil

#LOAD_LEVEL = 'LL' #Input low load (LL) or high load (HL)
# Go to line 1020 to see the file path.

def read_bus_number(case): #reading bus numbers
##    redirect.psse2py()
    psspy.psseinit()
##    psspy.psseinit(100000)
    psspy.case(case)
    ierr, iarray = psspy.abusint(-1, string='Number')
    Bus_number=iarray[0]
    return Bus_number



def read_bus_kv(case,bus): #reading bus numbers
##    redirect.psse2py()
    psspy.psseinit()
##    psspy.psseinit(100000)
    psspy.case(case)
    ierr, iarray = psspy.abusreal(-1, string='BASE')
    Bus_kv_tot=iarray[0]
    ierr, iarray = psspy.abusint(-1, string='Number')
    number=iarray[0]
    Bus_kv=Bus_kv_tot[number.index(bus)]
    return Bus_kv



#######################################################################################################################################
def read_bus_type(case): #reading bus types
    psspy.psseinit(100000)
    psspy.case(case)
    ierr, iarray = psspy.abusint(-1, string='Type')
    code=iarray[0]
    return code
#######################################################################################################################################
def read_machine_id(case): #reading bus types
    psspy.psseinit(100000)
    psspy.case(case)
    ierr, iarray = psspy.amachchar(-1, string='ID')
    ierr, iarray1 = psspy.amachint(-1,  string='NUMBER')
    machine_id=[x.strip(' ') for x in iarray[0]]
    machine_number=iarray1[0]
    machinedict=dict(zip(machine_number,machine_id))
    return machinedict


#######################################################################################################################################
def Branch_ID():  #read Branch ID
    ierr, br_id = psspy.abrnchar(-1,  string='ID')
    return br_id[0]
#######################################################################################################################################
def Tran_ID(): #read Transformer ID
    ierr, tran_id = psspy.atrnchar(-1,  string='ID')
    return tran_id[0]
#######################################################################################################################################
def Branch_Busnumbers():    #Find Branch buses
    ierr, branchfrom = psspy.abrnint(-1,  string='FROMNUMBER')
    ierr, branchto = psspy.abrnint(-1,  string='TONUMBER')
    return branchfrom[0], branchto[0]
#######################################################################################################################################
def Tran_Busnumbers():    #Find Tranfo buses
    ierr, tranfrom =psspy.atrnint(-1,  string='FROMNUMBER')
    ierr, tranto = psspy.atrnint(-1,  string='TONUMBER')

##    ierr, wind1 = psspy.atr3int(-1, string='WIND1NUMBER')
##    ierr, wind2 = psspy.atr3int(-1, string='WIND2NUMBER')
##    ierr, wind3 = psspy.atr3int(-1, string='WIND3NUMBER')

    return tranfrom[0], tranto[0]
#######################################################################################################################################
def Tran3_Busnumbers():    #Find 3winding Tranfo buses

    ierr, wind1 = psspy.atr3int(-1, string='WIND1NUMBER')
    ierr, wind2 = psspy.atr3int(-1, string='WIND2NUMBER')
    ierr, wind3 = psspy.atr3int(-1, string='WIND3NUMBER')


    return wind1[0], wind2[0], wind3[0]
#######################################################################################################################################
def isNaN(string):          #check the variable is None or not
    return string != string
#######################################################################################################################################
def dirCreateClean(path,fileTypes):        #creat a directory and remove existing files with specific file types
    def subCreateClean(subpath,fileTypes):
        try:
            os.mkdir(subpath)
        except OSError:
            pass

        #delete all files listed in fileTypes in the directory
        os.chdir(subpath)
        for type in fileTypes:
            filelist = glob.glob(type)
            for f in filelist:
                os.remove(f)
    currentDir = os.getcwd()
    subCreateClean(path,fileTypes)
##    subCreateClean(path+'\\Results',fileTypes)
##    subCreateClean(path+'\\PSSEOut',fileTypes)
    os.chdir(currentDir)
#######################################################################################################################################
#######################################################################################################################################
#%% Run Load flow

##def run_LoadFlow():
##    psspy.fnsl([0,0,0,1,1,0,0,0])
##    psspy.fnsl([0,0,0,1,1,0,0,0])
##    psspy.fnsl([0,0,0,1,1,0,0,0])
##    psspy.fnsl([0,0,0,1,1,0,0,0])
##    psspy.fnsl([2,0,0,1,1,0,0,0])
##    psspy.fnsl([2,0,0,1,1,0,0,0])
##    psspy.fdns([1,0,0,1,1,1,0,0])
##    psspy.fdns([1,0,0,1,1,1,0,0])
##    psspy.fdns([1,0,0,1,1,0,0,0])
##    psspy.fdns([1,0,0,1,1,0,0,0])
##    psspy.fnsl([1,0,0,1,1,0,0,0])	#Full Newton-Raphson
##    psspy.fnsl([1,0,0,1,1,0,0,0])
##    psspy.fnsl([1,0,0,1,1,0,0,0])
##    blownup = psspy.solved()
##	if blownup == 0:
##		return 0
##	else:
##		return 1

def run_LoadFlow():
##	psspy.fdns([1,0,0,1,1,1,0,0])
##	psspy.fdns([1,0,0,1,1,1,0,0])
##	psspy.fdns([1,0,0,1,1,0,0,0])
##	psspy.fdns([1,0,0,1,1,0,0,0])
    psspy.fnsl([0,0,0,1,1,0,0,0])
    psspy.fnsl([0,0,0,1,1,0,0,0])
    psspy.fnsl([0,0,0,1,1,0,0,0])
    psspy.fnsl([0,0,0,1,1,0,0,0])
    psspy.fnsl([2,0,0,1,1,0,0,0])
    psspy.fnsl([2,0,0,1,1,0,0,0])
    # psspy.fdns([1,0,0,1,1,1,0,0])
    # psspy.fdns([1,0,0,1,1,1,0,0])
    # psspy.fdns([1,0,0,1,1,0,0,0])
    # psspy.fdns([1,0,0,1,1,0,0,0])
    psspy.fnsl([1,0,0,1,1,0,0,0])	#Full Newton-Raphson
    psspy.fnsl([1,0,0,1,1,0,0,0])
    psspy.fnsl([1,0,0,1,1,0,0,0])
    blownup = psspy.solved()
    if blownup == 0:
        return 0
    else:
        return 1
#######################################################################################################################################
#######################################################################################################################################
#%% Prony Analysis

def prony(t, F, m):
	N= len(t)
	Amat = np.zeros((N-m, m))
	bmat = F[m:N]

	for jcol in range(m):
		Amat[:, jcol] = F[m-jcol-1:N-1-jcol]

	sol = np.linalg.lstsq(Amat, bmat)
	d = sol[0]

	# Solve the roots of the polynomial in step 2
	# first, form the polynomial coefficients
	c = np.zeros(m+1)
	c[m] = 1.
	for i in range(1,m+1):
		c[m-i] = -d[i-1]

	u = poly.polyroots(c)
	b_est = np.log(u)/(t[1] - t[0])

	# Set up LLS problem to find the "a"s in step 3
	Amat = np.zeros((N, m))
	bmat = F

	for irow in range(N):
		Amat[irow, :] = u**irow

	sol = np.linalg.lstsq(Amat, bmat)
	a_est = sol[0]

	return a_est, b_est
#######################################################################################################################################
#######################################################################################################################################



#%% Processing the signal
def processing_data(signal,time):
    status={}
    tol=0.01
    denOrder=50
    time=time.reset_index(drop=True)    #reset index of a signal
    signal=signal.reset_index(drop=True)
    initial=signal[0]
    Last=signal.iloc[-1]
    signal_stat=signal.describe()      #create statics of a signal
    if abs(Last-initial) < 0.01*signal_stat['mean']:
        status['error between last and begin point']=['ok']
    else:
        status['error between last and begin point']=[str(abs(Last-initial))]



    if signal_stat['std']<0.001:
        status['std']=['ok']


    if initial!=0 and Last==0:
        print('There is a trip')
        status['trip']=['should check']

    [Aset,Bset]=prony(time,signal,denOrder)
    undamp=0
    for k in range(len(Bset)):
        if Bset[k].real >0 and Bset[k].imag!=0:
            undamp+=1

    if undamp!=0:
        status['un-damped oscillation']=['should check']
    else:
        status['un-damped oscillation']=['ok']

    return status

##    Peak_loc_pos=scipy.signal.find_peaks(signal, height=signal_stat['mean']+signal_stat['mean']*tol)
##    Peak_loc_neg=scipy.signal.argrelextrema(np.array(signal_post),np.less)

##    scipy.signal.stft(signal)
#######################################################################################################################################
#######################################################################################################################################
#%% Find channels
def Find_Chan(BUS,QTY,Chan_Id):
    Ch,Ch_Id = [],[]

    for ich,ch in enumerate(BUS):
##        print ich,ch
        flag = 0
        Lastcheck=[]
        for id,chan in Chan_Id.items():

##            print id,chan
##            chan = "".join(chan.split())
            print (chan)

            if (str('VANGL') in str(QTY[ich])) :
                QTY_new='ANGL'
                check = re.findall(str(QTY_new).translate(str.maketrans('','', digits)), chan, flags=re.IGNORECASE)
            else:
                check = re.findall(str(QTY[ich]).translate(str.maketrans('','', digits)), chan.translate(str.maketrans('','', digits)), flags=re.IGNORECASE)
            busnumcheck = re.findall('[0-9]+',chan,flags=0) #find all numbers in the channel string.

            check1=[int(busno) for busno in busnumcheck if int(busno)==BUS[ich]] #check if busnumcheck has BUS number of intereste.
##            print check1
            if (len(re.findall("\d+", QTY[ich]))> 0) and  len(check1)>0:         #check ID in Channels
                if re.findall("\d+", QTY[ich])[0]==(busnumcheck[len(busnumcheck)-1]):
                    checkID=[1]
                else:
                    checkID=[]
            else:
                checkID=[1]

            if len(check)>0 and len(check1)>0 and chan.find(' TO ')<0 and len(checkID)>0: #chan.find('TO') to avoid finding channel for the power flow on the branch
               print (chan +'  '+str( id )+' '+ 'is matching for bus {}'.format(ch))
                ########
               print(f"digits = {digits}")
               if str('VANGL')==str(QTY[ich]):
                    if 'VOLT' in Chan_Id.items()[id-1][1]:

                       Ch.append(Chan_Id[id])
                       Ch_Id.append(id)
                       flag = 1
               elif str('ANGL')==str(QTY[ich]).translate(str.maketrans('','', digits)):

                    if 'VOLT' not in  Chan_Id.items()[id-1][1]:
                       Ch.append(Chan_Id[id])
                       Ch_Id.append(id)
                       flag = 1
               elif str('ANGL') not in str(QTY[ich]):
                    Ch.append(Chan_Id[id])
                    Ch_Id.append(id)
##               print Ch, Ch_Id
                    flag = 1

            elif not isinstance(BUS[ich],int):
                chan = "".join(chan.split())
##                print str('POWR'+str(ch[0])+'TO'+str(ch[1])+'CKT'+"'"+str(ch[2])+"'")
                if str('POWR'+str(ch[0])+'TO'+str(ch[1])+'CKT'+"'"+str(ch[2])+"'") in chan:
                    if (str(QTY[ich]) == 'P_FLOW'):
                       print (chan +'  '+str( id ) + 'is matching POWER FLOW for bus {}'.format(ch))
                       Ch.append(Chan_Id[id])
                       Ch_Id.append(id)
##                       print Ch, Ch_Id
                       flag = 1
                elif str('VARS'+str(ch[0])+'TO'+str(ch[1])+'CKT'+"'"+str(ch[2])+"'") in chan:

                    if (str(QTY[ich]) == 'Q_FLOW'):
                       print (chan +'  '+str( id ) + 'is matching VAR FLOW for bus {}'.format(ch))
                       Ch.append(Chan_Id[id])
                       Ch_Id.append(id)
##                       print Ch, Ch_Id
                       flag = 1
        if flag == 0:
            print ("Channel not found for bus {}".format(ch))
            Ch.append(None)
            Ch_Id.append(None)
##            print "Channel {} not found".format(ch)
##            print Ch,Ch_Id
##            print "*" * 10

##    print zip(Ch,Ch_Id)
    return Ch_Id,Ch

#######################################################################################################################################
#######################################################################################################################################

#%% Change Bus Number
def change_bus_number(sav_folder,df_Config_Data,OPDMS_connection,OPDMS_con_name):
    print("========================= CHANGING BUS NUMBERS =========================")
    Report={}
    df=pd.DataFrame.from_dict(Report)
    Number=df_Config_Data.loc[0,"Starting number"]
    Numbernew=int(Number)

    for caseid,case in enumerate(sav_folder):
        savecase=[]
        ad=[]
        tuples=[(str(case),'old Bus'),(str(case),'new Bus')]
        index = pd.MultiIndex.from_tuples(tuples)        #create index for each case bus number

        Action=str(df_Config_Data.loc[caseid,"Action"])

        Folder=CdirPath1+'\\'+str(case)+'\\'
        for address in os.listdir(Folder):                 #search in subfolders
            if os.path.isdir(os.path.join(Folder, address)) and address!='slackoff':
                for file in os.listdir(os.path.join(Folder, address)):
                    if (not file.endswith(".sav")):
                        continue
                    savecase.append(file)
                    ad.append(os.path.join(Folder, address))
        for file in os.listdir(Folder):
            if (not file.endswith(".sav")) or (file.endswith("_new.sav")):
                continue
            savecase.append(file)
            ad.append(Folder)
        if len(savecase)==1:
            Bus_old=read_bus_number(ad[0]+'\\'+savecase[0])
            Bus_types=read_bus_type(ad[0]+'\\'+savecase[0])
        else:
            raise Exception('There is no .sav case or more than one .sav cases in '+Folder)
        slack_old= Bus_old[Bus_types.index(3)]  #find slack bus
        Bus_new=[]
        counter=0
        if Action=='Yes':

            #########################################################################
            for k in range(len(Bus_old)):                                                     #Update bus numbers
                psspy.bus_chng_4(Bus_old[k],0,[_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f],_s)
                if len(eval(OPDMS_connection[caseid]))==1:
                    if Bus_old[k]==slack_old:
                        busnew_number=eval(OPDMS_connection[caseid])[0]
                        counter=1
                    else:
                        busnew_number=Numbernew-k+counter
                else:
                    busnew_number=Numbernew-k

                psspy.bus_number(Bus_old[k],busnew_number)
                Bus_new.append(busnew_number)
            Numbernew=Numbernew-k-1+counter
##            dirCreateClean(CdirPath1+'\\'+case+'\\',["*.sav","*.DAT"])

            branchfrom,branchto=Branch_Busnumbers()                                 #Change bus name
            Tranfrom,Tranto=Tran_Busnumbers()
            slack_new= Bus_new[Bus_types.index(3)]
            Gen_bus=[Bus_new[i] for i,val in enumerate(Bus_types) if val==2]
            if ((slack_new in (Tranfrom))) and (slack_new not in (branchfrom)) and (slack_new not in (branchto)) :
                impbus=Tranto[Tranfrom.index(slack_new)]
            elif (slack_new in (Tranto)) and (slack_new not in (branchfrom)) and (slack_new not in (branchto)) :
                impbus=Tranfrom[Tranto.index(slack_new)]
            elif(slack_new in (branchfrom)):
                impbus=branchto[branchfrom.index(slack_new)]

            elif(slack_new in (branchto)):
                impbus=branchfrom[branchto.index(slack_new)]

            bus_name=''
            txt=case.split('_')[1].split()
            if '&' in txt:
                txt.remove('&')
            if 'and' in txt:
                txt.remove('and')
            for k in range(min(len(txt),4)):
                if k==0:
                    bus_name = bus_name + txt[k][0]
                else:
                    bus_name = bus_name + ' ' + txt[k][0]
            
            psspy.bus_chng_4(impbus,0,[_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f],"".join(bus_name.split()))
            if not isNaN(OPDMS_con_name[caseid]):
                psspy.bus_chng_4(slack_new,0,[_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f],OPDMS_con_name[caseid])
            for k in range(len(Gen_bus)):
                psspy.bus_chng_4(Gen_bus[k],0,[_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f],"".join(bus_name.split())+' g'+str(k+1))
            newcase=savecase[0].replace(".sav",'')
            psspy.save(Folder+'\\'+newcase+'_new.sav')
        else:  #use old case number if new case already created just used that case number as old
                for file1 in os.listdir(Folder):
                    if (not file1.endswith("_new.sav")):
                        continue
                    newcase=file1
                if len(newcase)>0:
                    Bus_old=read_bus_number(Folder+'\\'+newcase)

                Bus_new=Bus_old



        ########################################################################
        df_new = pd.DataFrame([Bus_old,Bus_new],index=index)                     #update reprot dataframe
        df=df.append(df_new)

        ########################################################################
##

    return df


#######################################################################################################################################
#######################################################################################################################################






#%% Flat run
def Flat_Run(sav_folder,df_Config_Data,dynamic_sol_param):
    print("========================= RUNNING FLATRUN =========================")
        # Dynamic Solution Parameters --------------------------------------------------------
    
    max_solns = int(dynamic_sol_param.loc['max_solns','Var_Val'])                # network solution maximum number of iterations
    sfactor = float(dynamic_sol_param.loc['sfactor','Var_Val'])                    # acceleration sfactor used in the network solution (0.3 to 1)
    con_tol = float(dynamic_sol_param.loc['con_tol','Var_Val'])                # convergence tolerance used in the network solution (1usec)
    dT = float(dynamic_sol_param.loc['dT','Var_Val'])                       # simulation step time (1 to 4msec)
    frq_filter = float(dynamic_sol_param.loc['frq_filter','Var_Val'])              # filter time constant used in calculating bus frequancy deviations (4 to 16 msec, as long as this is at least 4 times time step)
    int_delta_thrsh = float(dynamic_sol_param.loc['int_delta_thrsh','Var_Val'])           # intermediate simulation mode time step threshold used in extended term simulations
    islnd_delta_thrsh = float(dynamic_sol_param.loc['islnd_delta_thrsh','Var_Val'])        # large (island frequency) simulation mode time step threshold used in extended term simulations
    islnd_sfactor = float(dynamic_sol_param.loc['islnd_sfactor','Var_Val'])              # large (island frequency) simulation mode acceleration factor used in extended term simulations
    islnd_con_tol = float(dynamic_sol_param.loc['islnd_con_tol','Var_Val'])           # large (island frequency) simulation mode convergence tolerance used in extended term simulations
    #Run and Plot Options
    AutoGNET=int(dynamic_sol_param.loc['AutoGNET','Var_Val'])                      #For missing machine models
    DisplayNonConvMon=int(dynamic_sol_param.loc['DisplayNonConvMon','Var_Val'])              #Display Network Convergerence Monitor
    PrintEvery=int(dynamic_sol_param.loc['PrintEvery','Var_Val'])                   #time steps
    WritetEvery=int(dynamic_sol_param.loc['WritetEvery','Var_Val'])                  #time steps
    PlotEvery=int(dynamic_sol_param.loc['Plt_t','Var_Val'])                     #time steps
    PSSE_case=[]
    initial=[]
    Report={}
    df=pd.DataFrame.from_dict(Report)
    for caseid,case in enumerate(sav_folder):
        Action=str(df_Config_Data.loc[caseid,"Action"])
        Gen_bus=[]
        df_new=[]
        if Action=='Yes':
            initial=[]
            PSSE_case=[]
            Folder=CdirPath1+'\\'+case+'\\'
            PSSE_case.append(case)
            #Input files
            for file in os.listdir(Folder):
                if (not file.endswith(".sav")):
                    continue
                savecase = file

            for file in os.listdir(Folder):
                if (not file.endswith(".dyr")):
                    continue
                dyrfile = file
            dll_models = []
            for file in os.listdir(Folder):
                if (not file.endswith(".dll")):
                    continue
                dll_models.append(file)
            System_bus=read_bus_number(Folder+savecase)
            Bus_types=read_bus_type(Folder+savecase)
            slack_bus= System_bus[Bus_types.index(3)]  #find slack bus
##            Gen_bus= [System_bus[Bus_types.index(2)] for x in ]
            Gen_bus=[System_bus[i] for i,val in enumerate(Bus_types) if val==2]                # find Generator buses
##            indexes = np.where(np.array(Bus_types) == 2)[0]
            print(f"((((((((((((((((caseid = {caseid}, case = {case}, folder = {Folder}")
            #Output File names
            CNVname = savecase.replace(".sav", "_CNV.cnv")
            SNPname = dyrfile.replace(".dyr", "_SNP.snp")
            OUTname = savecase.replace(".sav", "_Flatrun.out")
            PDEVname = savecase.replace(".sav", "_PDEV.txt")


            # Flat run config
            FlatRuntime = float(df_Config_Data.loc[caseid,"FlatRuntime(s)"])                #sec
            TotalSimulationDuration = float(df_Config_Data.loc[caseid,"TotalSimulationDuration(s)"])    #sec
            FaultDuration = float(df_Config_Data.loc[caseid,"FaultDuration(s)"])          #Duration of fault in sec

            # Disturbance test after Flat run

            FaultImpUnit = 1                #[1=MVA(Admittance), 2=Mho(Admittance), 3=Ohm(Impedance)]
            FaultImpedance = [0.0,-0.2E+10] #Admittance if selected 1 or 2 as units
            FaultBusNo=slack_bus


            FaultBuskV=read_bus_kv(Folder+savecase,slack_bus)
            psspy.close_powerflow()

            psspy.pssehalt_2()
            #PSS/E Initialization
            psspy.psseinit(100000)
            psspy.case(Folder+savecase)

            #Progress, Prompt, Alert, Report Settings
            psspy.lines_per_page_one_device(1,90)
            psspy.progress_output(2,Folder+PDEVname,[0,0])         #[2,0] to append on the existing file
            #psspy.alert_output(2,AlertFilename,[0,0])
            #psspy.prompt_output(2,ODEVname,[0,0])
            #psspy.report_output(2,ReportFilename,[0,0])

            #Setting up converted case for Dynamic Studies
            psspy.fnsl([1,0,0,1,1,0,0,0])
            psspy.cong(0)
            psspy.conl(0,1,1,[0,0],[ 100.0,0.0,0.0, 100.0])
            psspy.conl(0,1,2,[0,0],[ 100.0,0.0,0.0, 100.0])
            psspy.conl(0,1,3,[0,0],[ 100.0,0.0,0.0, 100.0])
            psspy.ordr(0)
            psspy.fact()
            psspy.fact()
            psspy.tysl(0)
            psspy.save(Folder+CNVname)

            #Load Dynamic Model file and dll libraries
            psspy.dyre_new([1,1,1,1],Folder+dyrfile,"","","")
            for dll in dll_models:
                psspy.addmodellibrary(Folder+dll)

            # Setting-up Dynamic Study Options and Solution Parameters
            psspy.set_netfrq(1)                 #Network Frequency Depedant
            psspy.set_relang(1,0,"")            #Set relative machine angles relative to system average angle
            psspy.dynamics_solution_param_2(intgar1=max_solns, realar1=sfactor, realar2 =con_tol, realar3=dT, realar4 =frq_filter, realar5=int_delta_thrsh, realar6=islnd_delta_thrsh, realar7=islnd_sfactor, realar8=islnd_con_tol)

            #Setting-up Channels for plotting
            psspy.change_channel_out_file(Folder+OUTname)
            psspy.chsb(0,1,[-1,-1,-1,1,1,0])    #M/C Angle
            psspy.chsb(0,1,[-1,-1,-1,1,2,0])    #M/C Pelec
            psspy.chsb(0,1,[-1,-1,-1,1,3,0])    #M/C Qelec
            psspy.chsb(0,1,[-1,-1,-1,1,4,0])    #M/C Eterm
            psspy.chsb(0,1,[-1,-1,-1,1,12,0])   #Bus Freq
            psspy.chsb(0,1,[-1,-1,-1,1,13,0])   #Bus Volt
            psspy.chsb(0,1,[-1,-1,-1,1,16,0])   #Branch FLow (P&Q)
            #Initialize and run upto Flatrun duration

            psspy.snap([-1,-1,-1,-1,-1],Folder+SNPname)
            psspy.case(Folder+CNVname)
            psspy.rstr(Folder+SNPname)
            psspy.strt_2([AutoGNET,DisplayNonConvMon],Folder+OUTname)
            psspy.strt(0,Folder+OUTname)
            psspy.strt(0,Folder+OUTname)
            psspy.strt(0,Folder+OUTname)
            psspy.run(DisplayNonConvMon,FlatRuntime,PrintEvery,WritetEvery,PlotEvery)
            initial_suspect = psspy.okstrt()
            #Apply Disturbance and run to end upto fault duration time
            psspy.dist_bus_fault(FaultBusNo,FaultImpUnit,FaultBuskV,FaultImpedance)
            psspy.change_channel_out_file(Folder+OUTname)
            psspy.run(DisplayNonConvMon,FlatRuntime+FaultDuration,PrintEvery,WritetEvery,PlotEvery)

            #Clear Disturbance and run upto end of simulation time
            psspy.dist_clear_fault(1)
            psspy.change_channel_out_file(Folder+OUTname)
            psspy.run(DisplayNonConvMon,TotalSimulationDuration,PrintEvery,WritetEvery,PlotEvery)
            if initial_suspect==0:
                initial.append('OK')
            else:
                initial.append('should check')
            # #save all progress output to pdv file
            # psspy.lines_per_page_one_device(2,60)
            # psspy.progress_output(1,"",[0,0])
            psspy.pssehalt_2()
        ########################################################################
            df_new = pd.DataFrame({'Case':PSSE_case})                     #update reprot dataframe
            df_new.insert(1,'Initial',initial)
        if len(df_new)!=0:
            df=df.append(df_new)

        ########################################################################
    return df
#######################################################################################################################################
#######################################################################################################################################
#%% Creating raw and sequence files

def raw_seq(sav_folder,df_Config_Data):
    print("========================= CREATING RAW AND SEQUENCE FILES =========================")


    for caseid,case in enumerate(sav_folder):
        lineinfo=df_Config_Data.loc[caseid,"System impedance Branch"]
        Action=str(df_Config_Data.loc[caseid,"Action"])
        zone=int(df_Config_Data.loc[caseid,"Zone"])
        Gen_bus=[]
        if Action=='Yes':
            Folder=CdirPath1+'\\'+case+'\\'

            #Input files
            for file in os.listdir(Folder):
                if (not file.endswith(".sav")):
                    continue
                savecase = file
                psspy.psseinit(100000)
                psspy.case(Folder+savecase)

            System_bus=read_bus_number(Folder+savecase)
            Bus_types=read_bus_type(Folder+savecase)
            slack_bus= System_bus[Bus_types.index(3)]  #find slack bus



         #########################################################################

            #remove slack bus gen
            machinedict=read_machine_id(Folder+savecase)
            psspy.machine_chng_2(slack_bus,str(machinedict[slack_bus]),[0,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])
            psspy.purgmac(slack_bus,str(machinedict[slack_bus]))
            psspy.bus_chng_4(slack_bus,0,[1,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f],_s)
         # Change Zone:
            for ib,bus in enumerate(System_bus):

                psspy.bus_chng_4(bus,0,[_i,int(zone),_i,_i],[_f,_f,_f,_f,_f,_f,_f],_s)
##                psspy.bus_chng_4(999791,0,[_i,3,_i,_i],[_f,_f,_f,_f,_f,_f,_f],_s)

        #########################################################################
        # remove line impedance
            if not isNaN(lineinfo):
                psspy.branch_data_3(eval(lineinfo)[0],eval(lineinfo)[1],str(eval(lineinfo)[2]), realar1=0,realar2=0.0001)
        #########################################################################
        # Remove title
            psspy.case_title_data("","")
            # create raw and seq files
            dirCreateClean(Folder+'\\'+'slackoff'+'\\',["*.sav","*.DAT"])
            newcase=savecase[0].replace(".sav",'')
            psspy.save(Folder+'\\'+'slackoff'+'\\'+newcase+'_infoff.sav')
            psspy.rawd_2(0,1,[1,1,1,0,0,0,1],0,Folder+savecase.replace(".sav",'')+'.raw')
            psspy.rwsq_2(0,1,[1,1,1,0,1],0,Folder+savecase.replace(".sav",'')+".seq")



#######################################################################################################################################
#######################################################################################################################################


#
#%% Analyse out
def Analyse_out(sav_folder,df_Config_Data):
    print("========================= ANALYSING OUT FILE =========================")


    Analyse_Bus_QTY=['VOLT','POWR','VARS']
    Analyse_Line_QTY=['P_FLOW','Q_FLOW']
    for caseid,case in enumerate(sav_folder):
        Action=str(df_Config_Data.loc[caseid,"Action"])
        Gen_bus=[]
        if Action=='Yes':
            Report={}
            df=pd.DataFrame.from_dict(Report)
            QTY=[]
            BUS=[]
            Folder=CdirPath1+'\\'+case+'\\'

            #Input files
            for file in os.listdir(Folder):
                if (not file.endswith(".out")):
                    continue
                Out_File = file
            for file in os.listdir(Folder):
                if (not file.endswith(".sav")):
                    continue
                savecase = file

            System_bus=read_bus_number(Folder+savecase)
            Bus_types=read_bus_type(Folder+savecase)
            slack_bus= System_bus[Bus_types.index(3)]  #find slack bus
##            Gen_bus= [System_bus[Bus_types.index(2)] for x in ]
            Gen_bus=[System_bus[i] for i,val in enumerate(Bus_types) if val==2]
            for Qty_ind, gen in enumerate(Gen_bus):

                for i,q in enumerate(Analyse_Bus_QTY):
                    QTY.append(q)
                    BUS.append(gen)

            QTY.append(Analyse_Bus_QTY[0])
            BUS.append(slack_bus)
            branchfrom,branchto=Branch_Busnumbers()                                 #Find Line to measure
            Tranfrom,Tranto=Tran_Busnumbers()

            if ((slack_bus in (Tranfrom))) and (slack_bus not in (branchfrom)) and (slack_bus not in (branchto)) :
                impbus=Tranto[Tranfrom.index(slack_bus)]

                line_id1=Tran_ID()
                line_id=line_id1[Tranfrom.index(slack_bus)].strip()
                Line=[slack_bus,impbus,str(line_id)]
            elif (slack_bus in (Tranto)) and (slack_bus not in (branchfrom)) and (slack_bus not in (branchto)) :
                impbus=Tranfrom[Tranto.index(slack_bus)]

                line_id1=Tran_ID()
                line_id=line_id1[Tranto.index(slack_bus)].strip()
                Line=[impbus,slack_bus,str(line_id)]
            elif(slack_bus in (branchfrom)):
                impbus=branchto[branchfrom.index(slack_bus)]

                line_id1=Branch_ID()
                line_id=line_id1[branchfrom.index(slack_bus)].strip()
                Line=[slack_bus,impbus,str(line_id)]
            elif(slack_bus in (branchto)):
                impbus=branchfrom[branchto.index(slack_bus)]

                line_id1=Branch_ID()
                line_id=line_id1[branchto.index(slack_bus)].strip()
                Line=[impbus,slack_bus,str(line_id)]
            for i,q in enumerate(Analyse_Line_QTY):
                QTY.append(q)
                BUS.append(Line)
            # Flat run config
            FlatRuntime = float(df_Config_Data.loc[caseid,"FlatRuntime(s)"])                #sec
            TotalSimulationDuration = float(df_Config_Data.loc[caseid,"TotalSimulationDuration(s)"])    #sec
            FaultDuration = float(df_Config_Data.loc[caseid,"FaultDuration(s)"])          #Duration of fault in sec


            ChnfObj = dt.CHNF(Folder+Out_File)
            Sht_ttl,Chan_Id,Chan_Data = ChnfObj.get_data()
            DF_Chan_data=pd.DataFrame(Chan_Data)
            Pre_data=DF_Chan_data[DF_Chan_data['time']<=FlatRuntime]
            During_data=DF_Chan_data[(DF_Chan_data['time']>FlatRuntime)&(DF_Chan_data['time']<=(FlatRuntime+FaultDuration))] #find pre, during and post contingency data
            Post_data=DF_Chan_data[(DF_Chan_data['time']>(FlatRuntime+FaultDuration+0.1))&(DF_Chan_data['time']<=TotalSimulationDuration)]
            Sim_Time = Chan_Data['time']
            start_time = Sim_Time[0]
            finish_time = Sim_Time[-1]
            chan_range = ChnfObj.get_range()
            channel_scale = ChnfObj.get_scale()

       #create index for each case bus number

            chans,Ch = Find_Chan(BUS,QTY,Chan_Id)
            no_chan = len(chans)
            for index, chan_id in enumerate(chans):
                if chan_id==None:
                    print ("*"*10)
                    print ("WARNING: Code cannot find a channel and at least one channel is None please investigate in "+ Folder)
                    warnings.warn("WARNING: Code cannot find a channel and at least one channel is None please investigatein in " + Folder)
                    print ("*"*10)
                if chan_id!=None:
                    tuples=[(str(Ch[index]),'Pre stat'),(str(Ch[index]),'Post stat')]
                    index = pd.MultiIndex.from_tuples(tuples)
                    pre_stat=processing_data(Pre_data[chan_id],Pre_data['time'])
                    post_stat=processing_data(Post_data[chan_id],Post_data['time'])
                    df_new = pd.DataFrame([pre_stat,post_stat],index=index)
                    df=df.append(df_new)
##                scipy.signal.find_peaks(Chan_Data[chans[0]])
##                scipy.signal.stft(Chan_Data[chans[0]])
##                statistics.pstdev(Chan_Data[chans[0]])
##                statistics.variance(Chan_Data[chans[0]])
##            Ch_Plt = Chan_Data[Figs_Plots[figure_key][subplot_key][k]][::]
            df.to_excel(writer, sheet_name = str(case))
            writer.save()



#######################################################################################################################################
#######################################################################################################################################
def Sorting_Key(string):
    if string.endswith(".sav"):
        num = int(string.split('_V')[1].split('.')[0])
    else:
        num = 0
    return int(num)

#%% Adding RUGs to OPDMS

def Adding_RUG(sav_folder,df_Config_Data):
    print(f"========================= ADDING RUG TO OPDMS CASE =========================")

    for caseid,case in enumerate(sav_folder):
        Action=str(df_Config_Data.loc[caseid,"Action"])
        OPDMS_case=df_Config_Data.loc[0,"OPDMS Case"]
        OPDMS_connection = (df_Config_Data["OPDMS connection Point"]).tolist()
        folderName = df_Config_Data.loc[caseid, "Folder name"]
        Gen_bus=[]
        if (OPDMS_case not in os.listdir(CdirPath1)) or (len(os.listdir(CdirPath1+'\\'+ OPDMS_case))==0):


            dirCreateClean(CdirPath1+'\\'+ OPDMS_case+'\\',[])
##            psspy.save(CdirPath1+'\\'+ OPDMS_case+'\\'+OPDMS_case_new+'_V0.sav') #for sajjad copy file there instead of ppsy.save
##            OPDMS_case_new=OPDMS_case_new+'_V0.sav'
##            shutil.move(OPDMS_case+'.sav',CdirPath1+'\\'+ OPDMS_case+'\\'+OPDMS_case+'_V0.sav')
            shutil.copy(CdirPath1+'\\'+OPDMS_case+'.sav',CdirPath1+'\\'+ OPDMS_case+'\\'+OPDMS_case+'.sav')
            OPDMS_case_new=OPDMS_case+'.sav'
            print(f"OPDMS_case_new (if) = {OPDMS_case_new} ")
        else:
            OPDMS_case_new=sorted(os.listdir(CdirPath1+'\\'+ OPDMS_case), key = Sorting_Key)[-1]   #select last version of the model in the OPDMS folder



        if Action=='Yes':
            Report={}
            System_path=[]
            psspy.psseinit(100000)
            psspy.close_powerflow()
            df=pd.DataFrame.from_dict(Report)
            Folder=CdirPath1+'\\'+case+'\\'
            for file in os.listdir(Folder):
                if (not file.endswith(".sav")):
                    continue
                savecase = file
            for file in os.listdir(Folder):
                if (not file.endswith(".raw")):
                    continue
                raw_File = file
            for file in os.listdir(Folder):
                if (not file.endswith(".seq")):
                    continue
                seq_file = file
            System_bus=read_bus_number(Folder+savecase)
            Bus_types=read_bus_type(Folder+savecase)
            slack_bus= System_bus[Bus_types.index(3)]  #find slack bus
##            Gen_bus= [System_bus[Bus_types.index(2)] for x in ]
            Gen_bus=[System_bus[i] for i,val in enumerate(Bus_types) if val==2]                # find Generator buses
            machinedict=read_machine_id(Folder+savecase)    #read machines id and bus number
            ###########################
            branchfrom,branchto=Branch_Busnumbers()                                 #Find branches
            Tranfrom,Tranto=Tran_Busnumbers()
            win1,win2,win3=Tran3_Busnumbers()
            all2wind=zip(Tranfrom,Tranto)
            all3wind=zip(win1,win2,win3)
            Tran3from=win1+win1
            Tran3to=win2+win3
            Linefrom=branchfrom+Tranfrom+Tran3from
            Lineto=branchto+Tranto+Tran3to
            Lines=zip(Linefrom,Lineto)

            ###########################
            # Create the graph of the system to find all the paths
            G = nx.Graph()
            G.add_nodes_from(System_bus)
            G.add_edges_from(Lines)
            End_nodes=[x for x in G.nodes() if G.degree(x)==1]
            if G.degree(slack_bus)==1:
                End_nodes.remove(slack_bus)
            for i,traget in enumerate(End_nodes):
                paths = nx.all_simple_paths(G, source=slack_bus, target=traget)
                for path in paths:
                    print(list(path))
                    System_path.append(list(path))
##            a=[item for item, count in collections.Counter(a).items() if count > 1]
            Connection_bus=list(dict.fromkeys([System_path[i][1] for i in range(len(System_path))])) #find buses connected to Slack

##            if len(eval(OPDMS_connection[caseid]))==1:
##                    if Bus_old[k]==slack_old:
##                        busnew_number=eval(OPDMS_connection[caseid])[0]
##                        counter=1
##                    else:
##                        busnew_number=Numbernew-k+counter
##            else:
##                    busnew_number=Numbernew-k




            # Load the OPDMS case
            ###
            psspy.close_powerflow()
            psspy.psseinit(100000)
            psspy.case(CdirPath1+'\\'+OPDMS_case+'\\'+OPDMS_case_new)
            err = run_LoadFlow()
            if err!=0:
                raise Exception('The OPDMS case of'+ OPDMS_case_new +' blown up before adding'+ case)
            else:
                opdms_ID = Sorting_Key(OPDMS_case_new)
                OPDMS_case_new=(OPDMS_case_new.replace('_V{}'.format(int(opdms_ID)), '_V{}'.format((int(opdms_ID)+1))))
                psspy.save(CdirPath1+'\\'+ OPDMS_case+'\\'+OPDMS_case_new)
            
            if (len(eval(OPDMS_connection[caseid])))!=1:
                ierr, Vsl = psspy.busdat(eval(OPDMS_connection[caseid])[0],'BASE')
                psspy.ltap(eval(OPDMS_connection[caseid])[0],eval(OPDMS_connection[caseid])[1],r"""1""", eval(OPDMS_connection[caseid])[2]*0.01,slack_bus,"", Vsl)
                if len(eval(OPDMS_connection[caseid]))>3:
                    psspy.ltap(eval(OPDMS_connection[caseid])[3],eval(OPDMS_connection[caseid])[4],r"""1""", eval(OPDMS_connection[caseid])[5]*0.01,slack_bus+5000,"", Vsl)
                    psspy.join(slack_bus, slack_bus+5000, 0)
            err = run_LoadFlow()
            ierr, Vbus = psspy.busdat(System_path[0][0],'PU')     #Point of connection bus voltage and angle before loading raw and seq
            ierr, angbus = psspy.busdat(System_path[0][0],'ANGLED')
            psspy.read(0,Folder+raw_File)
            psspy.resq(Folder+seq_file)
##            err = run_LoadFlow()

            for n, buses in enumerate(Connection_bus):
                psspy.bsysinit(1)
                psspy.bsyso(1,int(buses))
                psspy.dscn(int(buses))
                ierr, isl_buses = psspy.tree(1, 0)
                while isl_buses != 0:
                        ierr, isl_buses = psspy.tree(2, 1)
##                psspy.tree(2,1)
##                psspy.tree(2,1)

            psspy.bus_chng_4(System_path[0][0],0,[_i,_i,_i,_i],[_f, Vbus,angbus,_f,_f,_f,_f],_s) #update pont of connection angle and voltage
            err = run_LoadFlow()
            if err!=0:
                psspy.close_powerflow()
                raise Exception('The OPDMS case of'+ OPDMS_case_new +' blown up After adding'+ case)
            reconected_buses=[slack_bus]
            for i in range(len(System_path)):
                for k, grow_bus in enumerate((System_path[i])):
                    if grow_bus not in reconected_buses:
                        iteration=0
                        angle=[-30,30,0]
                        all2wind=zip(Tranfrom,Tranto)
                        output_2wind=list(filter( lambda x:grow_bus in x and System_path[i][k-1] in x , all2wind)) #find if the conncetion between System_path[i][k-1] and grow bus is two winding
                        while iteration<3:                                          # Give iteration for two winding transformer we need to check the angles and try +-30 degree
                        #####

                            ierr, Vbus = psspy.busdat(System_path[i][k-1],'PU')
                            ierr, angbus = psspy.busdat(System_path[i][k-1],'ANGLED')

                            reconected_buses.append(grow_bus)

    ##                        if ((grow_bus in Tran3from) or (grow_bus in Tran3to)) and ((System_path[i][k-1] in Tran3to) or (System_path[i][k-1] in Tran3from)): #check 3 winding buses and needs for turn on the third bus

    ##                            for wind3ind, val in enumerate (all3wind):

                                # output = [tupleitem
                                #    for tupleitem in all3wind
                                #    if ((grow_bus in tupleitem and System_path[i][k-1] in tupleitem)
                                #       )]

                            all3wind=zip(win1,win2,win3)
                            output_3wind=list(filter( lambda x:grow_bus in x and System_path[i][k-1] in x , all3wind))
                            new_bus=[elem for val in output_3wind for elem in list(val) if elem!= grow_bus and elem!= System_path[i][k-1]]
                                ###########################################
                            for item in new_bus:
                                psspy.bsysinit(1)
                                psspy.bsyso(1,item)
                                psspy.recn(item)
                                reconected_buses.append(item)
                                psspy.bus_chng_4(item,0,[_i,_i,_i,_i],[_f, Vbus,angbus,_f,_f,_f,_f],_s)
                            for b in range(len(output_3wind)):
                                psspy.three_wnd_imped_chng_4(output_3wind[b][0],output_3wind[b][1],output_3wind[b][2],r"""1""",[_i,_i,_i,_i,_i,_i,_i,_i,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,Vbus,angbus],_s,_s)
                            psspy.recn(grow_bus)
                            if len(output_2wind)==0 and len(output_3wind)==0:
                                iteration=2
                            psspy.bus_chng_4(grow_bus,0,[_i,_i,_i,_i],[_f, Vbus,angbus+angle[iteration],_f,_f,_f,_f],_s)
                            if grow_bus in Gen_bus:
                                psspy.machine_chng_2(grow_bus,str(machinedict[grow_bus]),[0,_i,_i,_i,_i,_i],[_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f,_f])
                            err = run_LoadFlow()


                            if err==0:
                                psspy.save(CdirPath1+'\\'+ OPDMS_case+'\\'+OPDMS_case_new)
                                iteration=3
                                
                            if err!=0:
                                iteration=iteration+1
                                if iteration==3:
                                    raise Exception('The OPDMS case of '+ OPDMS_case_new +' blown up in '+ case + " when the Bus {}".format(grow_bus)+' is on')
                                psspy.close_powerflow()
                                psspy.psseinit(100000)
                                psspy.case(CdirPath1+'\\'+OPDMS_case+'\\'+OPDMS_case_new)
                                reconected_buses.remove(grow_bus)
                                for item in new_bus:
                                    reconected_buses.remove(item)

            txt_file_path = CdirPath1+'\\'+OPDMS_case+'\\'+'ChangeLog_addingRUGs.txt'
            curr_dt = datetime.datetime.now()
            date_time = curr_dt.strftime("%Y/%m/%d %H:%M:%S")
            
            if not os.path.exists(txt_file_path):
                with open(txt_file_path, "w") as file:
                    file.write(f"{date_time}\n \t - RUG of {Folder} has been added to OPDMS sav case {OPDMS_case_new}\n")
                    file.write("\n")
            else:
                with open(txt_file_path, "a") as file:
                    file.write(f"{date_time}\n \t - RUG of {Folder} has been added to OPDMS sav case {OPDMS_case_new}\n")
                    file.write("\n")
                        # print(grow_bus)







#######################################################################################################################################
#######################################################################################################################################
#%% Main

if __name__ == '__main__':
 

    parser = argparse.ArgumentParser()
    parser.add_argument('-L', '--load-level', required=False)
    args = parser.parse_args()
    
    if args.load_level is None:
        load_level = "to be entered"
        bad_input = True
        while bad_input:
            load_level = input("Enter HL or LL: ")
            if not ((load_level == "HL") or (load_level == "LL")):
                print("Wrong input.")
            else:
                bad_input = False
    else:
        load_level = args.load_level
            
    #CdirPath1=os.getcwd()
    CdirPath1 = f"G:\Bungaban\PSSE_NEM_Models\PSSE Models\Bungaban_PSSE_NEM_Model_v0.1\{load_level}"
    config_file = f"OPDMS_Config_{load_level}.xlsx"
    report_path = os.path.join(CdirPath1,'Report.xlsx')

    if not os.path.exists(report_path):
    # Create an empty Excel file if it doesn't exist
        writer = pd.ExcelWriter(report_path, engine='openpyxl')
        # Create a dummy DataFrame and save it to create a visible sheet
        pd.DataFrame().to_excel(writer, index=False, header=False)
        writer.save()
        writer.close()
    writer = pd.ExcelWriter(CdirPath1+'\\'+'Report.xlsx', mode='a',on_sheet_exists="replace",engine = 'openpyxl')
    #%% Checklist
    Checklist= pd.read_excel(CdirPath1+'\\'+config_file, sheet_name="Config", engine='openpyxl')
    Checklist_Keys = Checklist.columns.tolist() #Column headings
    Checklist.set_index(Checklist_Keys[0],inplace=True); #Setting first column as index
    Checklist_test= list(Checklist.index.values)
    run_checklist=[]

    for k in Checklist_test:                                                             #check if the the function is selected
        action = str(Checklist.loc[k,"Action"])
        if action == "Yes":
            run_checklist.append(k)
    if any("change_bus numbers" in Functions for Functions in run_checklist):
        df_Config_Data = pd.read_excel(CdirPath1+'\\'+config_file, sheet_name="Change Number", engine='openpyxl')
        df_Config_Data.dropna(subset=["Folder name"], inplace=True)
        sav_cases = list(df_Config_Data["Folder name"])
        OPDMS_connection = (df_Config_Data["OPDMS connection Point"]).tolist()
        OPDMS_connection_name = (df_Config_Data["OPDMS connection Point Name"]).tolist()
        result_bus=change_bus_number(sav_cases,df_Config_Data,OPDMS_connection,OPDMS_connection_name)
        workBook = writer.book                                                              #to remove an exisitng Bus Number_report sheet, other wise it creates new sheets
##        try:
##            workBook.remove(workBook['Bus Number_report'])
##        except:
##            print("Worksheet does not exist")
##        finally:
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M')
        result_bus.to_excel(writer, sheet_name = 'Bus Number_report'+''+st)
        writer.save()

    if any("Flat Run" in Functions for Functions in run_checklist):
        df_Config_Data = pd.read_excel(CdirPath1+'\\'+config_file, sheet_name="Flat Run", engine='openpyxl')
        dynamic_sol_param=pd.read_excel(CdirPath1+'\\'+config_file,"Dynamic Solution Parameters",0, engine='openpyxl')
        dynamic_sol_param_Keys=dynamic_sol_param.columns.tolist(); #Column headings
        dynamic_sol_param.set_index(dynamic_sol_param_Keys[0],inplace=True); #Setting first column as index
        sav_cases = list(df_Config_Data["Folder name"])
        result_bus=Flat_Run(sav_cases,df_Config_Data,dynamic_sol_param)
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M')

        result_bus.to_excel(writer, sheet_name = 'Flat Run_report'+''+st)
        writer.save()
    if any("Create seq and raw" in Functions for Functions in run_checklist):
        df_Config_Data = pd.read_excel(CdirPath1+'\\'+config_file, sheet_name="Create seq and raw", engine='openpyxl')
        df_Config_Data.dropna(subset=["Folder name"], inplace=True)
        sav_cases = list(df_Config_Data["Folder name"])
        raw_seq(sav_cases,df_Config_Data)
    if any("Analyse out file" in Functions for Functions in run_checklist):
        df_Config_Data = pd.read_excel(CdirPath1+'\\'+config_file, sheet_name="Flat Run", engine='openpyxl')
        sav_cases = list(df_Config_Data["Folder name"])
        result_analys=Analyse_out(sav_cases,df_Config_Data)
        ts = time.time()
        st = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M')

        result_analys.to_excel(writer, sheet_name = 'Analysis_report'+''+st)
        writer.save()

    if any("Add to OPDMS" in Functions for Functions in run_checklist):
        df_Config_Data = pd.read_excel(CdirPath1+'\\'+config_file, sheet_name="Add to OPDMS", engine='openpyxl')
        sav_cases = list(df_Config_Data["Folder name"])
        Adding_RUG(sav_cases,df_Config_Data)
##        ts = time.time()
##        st = datetime.datetime.fromtimestamp(ts).strftime('%Y%m%d_%H%M')

##        result_analys.to_excel(writer, sheet_name = 'Analysis_report'+''+st)
##        writer.save()
    writer.close()
##        print(caseid)
##        print(case)





