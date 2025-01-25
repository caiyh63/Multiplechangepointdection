# Multiple Change Point Dection
* The code package mainly apply to repeated the results conducted by Zhao et al. (2018). See Fig.2  
* Preparation: Before run the RJMCMC code, you should the prepare the file invoving the the time series of tropical cyclone genesis frequency over the North Atlantic  
  during 1979-2014 (ALTC.txt). The time series in file is taken as the optional input in built-in function file "MultipleChangePointDetection.m".
* Run the Matlab code "bayc_TC_test.m" that is  used to generate the Posterior probility of each candidate hypothesis (ProbHypothesis.txt) and sample of changepoint (SampleOfChangePoint1.txt). Note that it should need  built-in function file "MultipleChangePointDetection.m"  saved to your working directory.     
* The algorithm from "MultipleChangePointDetection.m" is developed in Zhao and Chu (2009)    
* Plot your results by running the python code "baychangep.py". Files "ALTC.txt", "ProbHypothesis.txt" and "SampleOfChangePoint1.txt" are used to read and plot as below:     
>>> a) Time Series of  NATL TCGF  
>>> b) Posterior probability of each candidate hypothesis  
>>> c) Posterior Probability Mass Funtion (PMF)  
* See "changepoint.jpg"
## Reference:
** Zhao H, Duan X, Raga GB, Sun F 2018: Potential large-scale forcing mechanisms driving enhanced North Atlantic tropical cyclone activity since the mid-1990s. J. Climate  
** Zhao, X. and P.S. Chu, 2009: Bayesian Change-Point Analysis for Extreme Events (Typhoons, Heavy Rainfall, and Heat Waves): A RJMCMC Approach, J. Climate  
 
