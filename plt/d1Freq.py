import numpy as np
import glob
import sys
import re
import matplotlib.pyplot as plt
from datetime import datetime
import json
import pandas as pd
import scipy.stats as stats
from decimal import Decimal, getcontext



#this script computes the distribution of the intracell distance between xAj and xBj


if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit()
rowNum=0#int(sys.argv[1])
unitCellNum=int(sys.argv[1])

csvDataFolderRoot="../dataAll/dataAllUnitCell"+str(unitCellNum)+"/row"+str(rowNum)+"/csvOutAll/"

inParamFileName="../V_1d_quartic_2min.csv"



inDf=pd.read_csv(inParamFileName)
oneRow=inDf.iloc[rowNum,:]
a=float(oneRow.loc["a"])
b=float(oneRow.loc["b"])
c=float(oneRow.loc["c"])
f=float(oneRow.loc["f"])

d1=float(oneRow.loc["d1"])
d2=float(oneRow.loc["d2"])

TVals=[]
TFileNames=[]

for TFile in glob.glob(csvDataFolderRoot+"/T*"):

    matchT=re.search(r"T(\d+(\.\d+)?)",TFile)
    # if float(matchT.group(1))<1:
    #     continue

    if matchT:
        TFileNames.append(TFile)
        TVals.append(float(matchT.group(1)))


sortedInds=np.argsort(TVals)
sortedTVals=[TVals[ind] for ind in sortedInds]
sortedTFiles=[TFileNames[ind] for ind in sortedInds]

def V1(x):

    val=(x-d1)**4-a*(x-d1)**2+b
    val*=c
    return val
def format_using_decimal(value, precision=4):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)
def plt_d1_hist(oneTFile):
    """

    :param oneTFile: corresponds to one temperature
    :return:
    """

    matchT=re.search(r'T([-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?)',oneTFile)
    TVal=float(matchT.group(1))
    TStr=format(TVal)

    U_distPath=oneTFile+"/U_dist/U_distData.csv"

    df=pd.read_csv(U_distPath)

    lattice_arr=np.array(df.iloc[:,1:])

    _,nCol=lattice_arr.shape

    N=int(nCol/2)

    xA_arr=lattice_arr[:,0:N]
    xB_arr=lattice_arr[:,N:]

    d1_arr=xB_arr-xA_arr

    d1_vec=d1_arr.reshape(-1)

    d1_left=d1-np.sqrt(a/2)*1.5

    d1_right=d1+np.sqrt(a/2)*1.5

    xValsAll=np.linspace(d1_left,d1_right,100)
    V1ValsAll=np.array([V1(x) for x in xValsAll])

    V1Mean=np.mean(V1ValsAll)
    # V1ValsAll-=V1Mean


    fig, ax1 = plt.subplots()

    ax1.hist(d1_vec, bins=100, density=True, alpha=0.6, color='red',label="occurrence")
    ax1.set_xlabel('Intracell distance')
    ax1.set_ylabel('Density', color='black')

    ax2 = ax1.twinx()
    ax2.plot(xValsAll, V1ValsAll, 'b-', lw=1,label="V1")
    ax2.set_ylabel("V1")

    plt.legend(loc="best")
    d1_HistOut="T"+str(TVal)+"d1_hist.png"
    plt.title("T="+TStr+", histogram of intracell distance")
    plt.savefig(oneTFile+"/"+d1_HistOut)

    plt.close()









for k in range(0,len(sortedTFiles)):
    print("processing T="+str(sortedTVals[k]))
    oneTFile=sortedTFiles[k]
    plt_d1_hist(oneTFile)