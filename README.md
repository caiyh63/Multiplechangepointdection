# Multiplechangepointdection
* The code package mainly apply to repeated the results conducted by Zhao et al. (2018). See Fig.2  

* "bayc_TC_test.m" is Matlab code used to wirte the Posterior probility of each candidate hypothesis (ProbHypothesis.txt) and sample of changepoint (SampleOfChangePoint1.txt) it should need to built-in function file "MultipleChangePointDetection.m", which save it to your working directory.   
# The algorithm from "MultipleChangePointDetection.m" is developed in Zhao and Chu (2009)
# The input file in bayc_TC_test.m is the time series of tropical cyclone genesis frequency over the North Atlantic during 1979-2014 (ALTC.txt), writted by NCL code (NATL.ncl)
# Files "ALTC.txt", "ProbHypothesis.txt" and "SampleOfChangePoint1.txt" are used to read in "baychangep.py" and plot as below: 
# a) Time Series of  NATL TCGF
# b) Posterior probability of each candidate hypothesis
# c) Posterior Probability Mass Funtion (PMF)
# See "changepoint.jpg"
# Reference:
# Zhao H, Duan X, Raga GB, Sun F 2018: Potential large-scale forcing mechanisms driving enhanced North Atlantic tropical cyclone activity since the mid-1990s. J. Climate
# Zhao, X. and P.S. Chu, 2009: Bayesian Change-Point Analysis for Extreme Events (Typhoons, Heavy Rainfall, and Heat Waves): A RJMCMC Approach, J. Climate
 
