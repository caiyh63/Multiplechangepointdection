# This code repeat the result from Zhao et al. (2018) (See Fig. 2) 
# that detect the change point of tropical cyclone genesis frequency over the North Atlantic (NATL TCGF)
# This deteced change point base on the Bayesian changepoint method (Zhao and Chu 2010)
# a) Time Series of  NATL TCGF
# b) Posterior probability of each candidate hypothesis
# c) Posterior Probability Mass Funtion (PMF)
# 
#
# Last Edited:
#    2023-10        First Create
# 
#  Yuhao CAI From SYSU

import numpy as np
import pandas as pd

ALTC    = np.loadtxt("./ALTC.txt")
Hprob   = np.loadtxt("./ProbHypothesis.txt")
SOCP    = np.loadtxt("./SampleOfChangePoint1.txt")

SOCP    = SOCP + 1979 -1
year    = np.arange(1979,2015,1)
freq    = np.zeros(len(year))

# Obtain the frequency and probality distribution
for yy in range(len(year)):
    freq[yy] = np.sum(SOCP == year[yy])
Prob    = freq/len(SOCP)
max_idx = np.where(Prob == max(Prob)) # Attain the index of maximum values
dd      = max_idx[0][0]
# average value during two sub-period
avgpr,avgpo   = np.zeros(dd), np.zeros(len(year)-dd)
avgpr[:]  = np.mean(ALTC[:dd])
avgpo[:]  = np.mean(ALTC[dd:])

# POLT
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator #从pyplot导入MultipleLocator类，这个类用于设置刻度间隔


fig   = plt.figure(figsize = (18, 12))
ax1   = fig.add_subplot(311)
ax2   = fig.add_subplot(312)
ax3   = fig.add_subplot(313)
ax1.set_position([0.3, 0.70, 0.3, 0.25])
ax2.set_position([0.3, 0.38, 0.3, 0.25])
ax3.set_position([0.3, 0.05, 0.3, 0.25])

#公共绘图资源设置--=======================================================
def line_plot(figs,strs):
    #设置坐标轴
    figs.set_xticks(year)
    figs.tick_params(labelsize=18)
    #设置坐标轴刻度
    figs.yaxis.set_major_locator(MultipleLocator(10)) #set interval of axis
    figs.yaxis.set_minor_locator(MultipleLocator(2))
    figs.xaxis.set_major_locator(MultipleLocator(5)) #set interval of axis
    figs.xaxis.set_minor_locator(MultipleLocator(1))
    figs.tick_params(axis='y',which='minor',direction='in',length=5,width=2) #调整刻度线长度
    figs.tick_params(axis='y',which='major',direction='in',length=5,width=2,pad=20) #调整刻度线长度
    figs.tick_params(axis='x',which='major',direction='in',length=5) #调整刻度线长度
    figs.spines['right'].set_linewidth(2)
    figs.spines['left'].set_linewidth(2)
    figs.spines['top'].set_linewidth(2)
    figs.spines['bottom'].set_linewidth(2)
    # 设置标题
    figs.set_title(strs,loc='left',fontsize =20)
    
# 1st subplot
ax1.set_ylim(0,25)
ax1.set_xlim(1979,2014)
line_plot(ax1,"a) Peak season NATL TC")
l1     = ax1.plot(year,ALTC,color='k',linewidth=2,marker='.',markersize=15)


l21 = ax1.plot(year[:dd],avgpr,color='grey',linewidth=2,ls='--')
l22 = ax1.plot(year[dd:],avgpo,color='grey',linewidth=2,ls='--')

x = range(7)

def bar_plot(ax, Lstrs):
  # set axis source
  ax.yaxis.set_major_locator(MultipleLocator(0.1)) #set interval of axis
  ax.tick_params(axis='y',which='minor',direction='in',length=5,width=2) #调整刻度线长度
  ax.tick_params(axis='y',which='major',direction='in',length=5,width=2,pad=20) #调整刻度线长度
  ax.tick_params(axis='x',which='major',direction='in',length=5) #调整刻度线长度
  ax.spines['bottom'].set_linewidth(2)            ###设置底部坐标轴的粗细
  ax.spines['left'].set_linewidth(2)              ####设置左边坐标轴的粗细
  ax.spines['right'].set_linewidth(2)             ###设置右边坐标轴的粗细
  ax.spines['top'].set_linewidth(2)               ####设置上部坐标轴的粗细
  ax.set_title(Lstrs,loc='left',fontsize=20,fontweight='light')
  ax.tick_params(labelsize=18)

# 2nd subplot
bar_plot(ax2, "b) Posterior Probability ")
ax2.set_ylim(0,0.71,0.1)
ax2.set_xlim(-1,7)    
ax2.set_xticks(x) # 设置刻度
ax2.set_xlabel("Number of Change Points",fontsize=15)
ax2.bar(x,Hprob,0.3,color='k',align='center',alpha=1.0)

# 3rd subplot
bar_plot(ax3, "c) Change Point")
ax3.set_ylabel("Posterior PMF",fontsize=15)
ax3.bar(year,Prob,0.8,color=np.where(Prob == max(Prob),'r','k'),align='center',alpha=0.7) # +bar_widh/4
ax3.text(year[max_idx]-1.5,Prob[max_idx]-0.01,str(year[max_idx[0][0]]),c='k',fontsize=15)


jpg_fig = plt.gcf() # 'get current figure
jpg_fig.savefig('changepoint.jpg', format='jpg', dpi=300)
plt.show()
