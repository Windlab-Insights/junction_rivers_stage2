import sys
import psse34
import psspy
from psspy import _i, _f, _s
import redirect
import pandas as pd

psspy.psseinit()
redirect.psse2py()
psspy.case(r"G:\Bungaban\PSSE_NEM_Models\PSSE Models\Bungaban_PSSE_NEM_Model_v0.1\HL\OPDMS Files\SpringHi-20220905-190113-SystemNormal.sav")
# psspy.case(r"G:\Bungaban\PSSE_NEM_Models\PSSE Models\Bungaban_PSSE_NEM_Model_v0.1\HL\NEMSavCases\ManualGeneratorDispatch\V07\ManualGeneratorDispatch_V0.7.sav")
# Get network data
ierr, BusArray0 = psspy.abusint(-1,2,"NUMBER")
ierr, BusArray1 = psspy.abuschar(-1, 2, "NAME")
ierr, BusArray2 = psspy.abusreal(-1,2,['BASE', "PU"])

ierr, BranchArray0 = psspy.abrnchar(-1, 2, 3, 2, 1,["ID",'FROMNAME','TONAME'])
ierr, BranchArray1 = psspy.abrnint(-1, 2, 3, 2, 1, ["FROMNUMBER",'TONUMBER'])
ierr, BranchArray2 = psspy.abrnreal(-1, 2, 3, 2, 1, ['PCTRATEA','PCTRATEB'])

ierr, TranArray0 = psspy.atrnchar(-1, 2, 3, 2, 1,["ID",'FROMNAME','TONAME'])
ierr, TranArray1 = psspy.atrnint(-1, 2, 3, 2, 1, ["FROMNUMBER",'TONUMBER'])
ierr, TranArray2 = psspy.atrnreal(-1, 2, 3, 2, 1, ['PCTRATEA','PCTRATEB','NOMV1','NOMV2'])

ierr, Tran3Array0 = psspy.awndchar(-1, 2, 3, 2, 1,["ID",'WIND1NAME','WIND2NAME','WIND3NAME'])
ierr, Tran3Array1 = psspy.awndint(-1, 2, 3, 2, 1, ["WIND1NUMBER",'WIND2NUMBER','WIND3NUMBER'])
ierr, Tran3Array2 = psspy.awndreal(-1, 2, 3, 2, 1, ['PCTRATEA','PCTRATEB'])

ierr, qni1 = psspy.brnflo(231694,460391,"1")
if qni1.real>=0:
    qni1Dir = "NSW->QLD"
else:
    qni1Dir = "QLD->NSW"
qni1 = abs(qni1)

ierr, qni2 = psspy.brnflo(231692,460392,"2")
if qni2.real>=0:
    qni2Dir = "NSW->QLD"
else:
    qni2Dir = "QLD->NSW"
qni2 = abs(qni2)

ierr, mud1 = psspy.brnflo(276440,440440,"1")
if mud1.real>=0:
    mud1Dir = "NSW->QLD"
else:
    mud1Dir = "QLD->NSW"
mud1 = abs(mud1)

ierr, mud2 = psspy.brnflo(276440,440440,"2")
if mud2.real>=0:
    mud2Dir = "NSW->QLD"
else:
    mud2Dir = "QLD->NSW"
mud2 = abs(mud2)

busData = pd.DataFrame({
    'Bus Number': BusArray0[0],
    'Bus Name': BusArray1[0],
    'Bus V_Base': BusArray2[0],
    'Bus V_PU': BusArray2[1]
})

branchData = pd.DataFrame({
    "CircuitID": BranchArray0[0],
    "FromName": BranchArray0[1],
    "ToName": BranchArray0[2],
    "FromNum": BranchArray1[0],
    "ToNum": BranchArray1[1],
    "RateALoading": BranchArray2[0],
    "RateBLoading": BranchArray2[1],
})  

tranData = pd.DataFrame({
    "CircuitID": TranArray0[0],
    "FromName": TranArray0[1],
    "ToName": TranArray0[2],
    "FromNum": TranArray1[0],
    "ToNum": TranArray1[1],
    "RateALoading": TranArray2[0],
    "RateBLoading": TranArray2[1],
    "VWIND1": TranArray2[2],
    "VWIND2": TranArray2[3],
})  

tran3Data = pd.DataFrame({
    "CircuitID": Tran3Array0[0],
    "WIND1NAME": Tran3Array0[1],
    "WIND2NAME": Tran3Array0[2],
    "WIND3NAME": Tran3Array0[3],
    "WIND1NUM": Tran3Array1[0],
    "WIND2NUM": Tran3Array1[1],
    "WIND3NUM": Tran3Array1[2],
    "RateALoading": Tran3Array2[0],
    "RateBLoading": Tran3Array2[1],
})


filt_busData = busData[busData['Bus V_Base'] >= 110]
sort_df = filt_busData.sort_values(by='Bus V_PU', ascending=False)
top_10 = sort_df.head(10)
low_10 = sort_df.tail(10).sort_values(by='Bus V_PU', ascending=True)
print("~~~~~~~~Top Ten Voltage Mag~~~~~~~~")
print(top_10[['Bus Name', 'Bus V_PU']])

print("~~~~~~~~Low Ten Voltage Mag~~~~~~~~")
print(low_10[['Bus Name', 'Bus V_PU']])

print("~~~~~~~~Top Ten Line Loading~~~~~~~")
sort_df = branchData.sort_values(by='RateALoading', ascending=False)
top_10 = sort_df.head(10)
print(top_10[['FromName', 'ToName',"RateALoading" ]])

print("~~~~~~~~Top Ten 2TX Loading~~~~~~~~")
sort_df = tranData.sort_values(by='RateALoading', ascending=False)
top_10 = sort_df.head(10)
print(top_10[['FromName', 'ToName',"RateALoading" ]])


print("~~~~~~~~Top Ten 3TX Loading~~~~~~~~")
sort_df = tran3Data.sort_values(by='RateALoading', ascending=False)
top_10 = sort_df.head(10)
print(top_10[['WIND1NAME',"WIND2NAME", "WIND3NAME","RateALoading" ]])


print("~~~~~~~~Interconnector Flow~~~~~~~~")
print(f"QNI1: {round(qni1,2)} (58.62), Dir: {qni1Dir} (NSW->QLD)")
print(f"QNI2: {round(qni2,2)} (58.62), Dir: {qni2Dir} (NSW->QLD)")
print(f"MUD1: {round(mud1,2)} (48.72), Dir: {mud1Dir} (QLD->NSW)")
print(f"MUD2: {round(mud2,2)} (48.72), Dir: {mud2Dir} (QLD->NSW)")

print("Stop")