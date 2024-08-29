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


#this script computes correlation function between d1 and d2


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



def cor_d1_d2(oneTFile,ind1,ind2):
    """

    :param oneTFile: corresponds to one temperature
    :param ind1: index of d1 in vec
    :param ind2: index of d2 in vec
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

    d1_vec_tmp=xB_arr[:,ind1]-xA_arr[:,ind1]

    d2_vec_tmp=xA_arr[:,ind2+1]-xB_arr[:,ind2]

    covTmp=np.cov(d1_vec_tmp,d2_vec_tmp,ddof=1)
    return covTmp

cov_d1d1_vec=[]
cov_d1d2_vec=[]
cov_d2d2_vec=[]

ind1=1
ind2=2
for k in range(0,len(sortedTFiles)):
    print("processing T="+str(sortedTVals[k]))
    oneTFile=sortedTFiles[k]

    covTmp=cor_d1_d2(oneTFile,ind1,ind2)
    cov_d1d1_tmp=covTmp[0,0]
    cov_d1d2_tmp=covTmp[0,1]
    cov_d2d2_tmp=covTmp[1,1]

    cov_d1d1_vec.append(cov_d1d1_tmp)
    cov_d1d2_vec.append(cov_d1d2_tmp)
    cov_d2d2_vec.append(cov_d2d2_tmp)

cov_d1d1_vec=np.array(cov_d1d1_vec)
cov_d1d2_vec=np.array(cov_d1d2_vec)
cov_d2d2_vec=np.array(cov_d2d2_vec)

sortedTVals=np.array(sortedTVals)
TInds=np.where(sortedTVals<100)
TToPlt=sortedTVals[TInds]


plt.figure()
plt.plot(TToPlt,cov_d1d1_vec[TInds],color="black",marker="o",lw=1)
plt.xlabel("$T$")

plt.ylabel("cov(d1,d1)")
plt.title("cov(d1,d1)")
plt.savefig(csvDataFolderRoot+"/d1d1_"+str(ind1)+"_"+str(ind2)+".png")
plt.close()

plt.figure()
plt.plot(TToPlt,cov_d1d2_vec[TInds],color="red",marker="x",lw=1)
plt.xlabel("$T$")
plt.ylabel("cov(d1,d2)")
plt.title("cov(d1,d2)")
plt.savefig(csvDataFolderRoot+"/d1d2_"+str(ind1)+"_"+str(ind2)+".png")
plt.close()


plt.figure()
plt.plot(TToPlt,cov_d2d2_vec[TInds],color="blue",marker=".",lw=1)
plt.xlabel("$T$")
plt.ylabel("cov(d2,d2)")
plt.title("cov(d2,d2)")
plt.savefig(csvDataFolderRoot+"/d2d2_"+str(ind1)+"_"+str(ind2)+".png")
plt.close()