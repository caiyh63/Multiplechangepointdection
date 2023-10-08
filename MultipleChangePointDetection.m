function [ProbHypothesis, SampleOfHypothesisIndex, SampleOfChangePoint, SampleOfRate] = MultipleChangePointDetection(Series, TimeIndex, Length_Burnin, Length_Main, NumberOfHypotheses, DisplayMode)

%MULTIPLECHANGEPOINTDETECTION Implement an RJMCMC algorithm to identify and quantify 
%   the multiple abrupt shifts of a given rare event count series through Bayesian 
%   Analysis. By default, this function provides the posterior probability of possible 
%   change-points and draw the resulted empirically fitted PMF bar chart; also it 
%   provides the estimated posterior PMF of the change-points under the winning hypothesis.
%
%   [ProbHypothesis, SampleOfHypothesisIndex, SampleOfChangePoint, SampleOfRate] =
%       MultipleChangePointDetection(Series, TimeIndex, Length_Burnin, Length_Main, NumberOfHypotheses, DisplayMode)
%   [ProbHypothesis] = MultipleChangePointDetection(Series)
%
%   Optional Inputs: TimeIndex, Length_Burnin, Length_Main, NumberOfHypotheses, DisplayMode
%
% Inputs:
%   Series - Vector of observations of a rare event count series.
%
% Optional Inputs:
%   TimeIndex - The time index vector of the observation. If it's not null and its length is not equal to the length 
%       of input count series, "Series", this function shall choose the first value of this vector as the starting index.
%       By default, when it's set as null, this function use [1:N] as series index, where N is the length of input series.
%
%   Length_Burnin - The length of the burn-in period. By default, it's set as 2000.
%
%   Length_Main - The length of the main period. By default, it's set as 10000.
%
%   NumberOfHypotheses - Number of all the candidate hypotheses. By default, it's
%       set as the minimum value between 10 and length of series divided by 4 (floor).
%
%   DisplayMode - By default, it's set as 1 and during the simulation, it will display
%       the number of the finished MCMC iterations in every 1000 counts. Also, in the 
%       end of the simulation, the function will draw the bar chart for the PMF for
%       each candidate hypothesis. If this value is 0, the function will show nothing. 
%
% Outputs:
%   ProbHypothesis - Estimated posterior probability of each of the candidate hypothesis.
%
%   SampleOfHypothesisIndex - The samples of accepted hypothesis index drawn from the main
%       period of the MCMC simulation.
%
%   SampleOfChangePoint - The samples of the change-points under the accepted hypothesis index 
%       drawn from the main period of the MCMC simulation. It's a 3-dimensional matrix in the
%       order ("Samples", "index of iterations", "hypothesis index"). The "samples" within 
%       each iteration is listed in an ascending order with respect to time, like (tao1, tao2,...). 
%       In each MCMC iteration, we pad zero if there is no sample applicable to a given hypothesis 
%       and/or given sample index. 
%
%   SampleOfRate - The samples of the rates under the accepted hypothesis index drawn from 
%       the main period of the MCMC simulation. It's also a 3-dimensional matrix in the
%       order ("Samples", "index of iterations", "hypothesis index"). The "samples" within 
%       each iteration is listed in an ascending order with respect to time, like (lbd1, lbd2,...). 
%       In each MCMC iteration, we pad zero if there is no sample applicable to a given hypothesis 
%       and/or given sample index.
%
% Examples:
%   [ProbHypothesis] = MultipleChangePointDetection(Series)
%
%   [ProbHypothesis] = MultipleChangePointDetection(Series, TimeIndex)
%
%   [ProbHypothesis, SampleOfHypothesisIndex, SampleOfChangePoint, SampleOfRate] =
%       MultipleChangePointDetection(Series, TimeIndex(1), 10000, 40000, 8, 0)
%
%   [ProbHypothesis, SampleOfHypothesisIndex, SampleOfChangePoint, SampleOfRate] =
%       MultipleChangePointDetection(Series, [], [], [], 10)
%
%   [ProbHypothesis, SampleOfHypothesisIndex] = MultipleChangePointDetection(Series, [], [], 4000)


%   Made by X. Zhao on 01-01-09
%   Reference: Zhao, X. and P.S. Chu, 2009: Bayesian Change-Point Analysis for Extreme Events 
%       (Typhoons, Heavy Rainfall, and Heat Waves): A RJMCMC Approach, J. Climate


%-----------------
% initialization 
%-----------------
if (nargin == 0)
    return;
end

% Reading Data
Y = (Series(:))';                               % the tested time series
N = length(Y);                                  % length of time series
% Setting System Parameters
timeIndex_default = [1:N];
buildin_default = 2000;                         % default length of burn-in period 
all_default = 10000;                            % default length of main period
MK_default = min(10, floor(N / 5));             % default total number of candidate hypotheses
DPMode_default = 1;                             % default display mode. 1: Display every 100 iterations; 0 : never display until the finish.
if (nargin >= 2) & ~isempty(TimeIndex)
    timeIndex = TimeIndex;
    if (length(TimeIndex) ~= N)
        timeIndex = TimeIndex(1) - 1 + [1:N];
    end
else
    timeIndex = timeIndex_default;
end
if (nargin >= 3) & ~isempty(Length_Burnin)
    buildin = Length_Burnin;
else
    buildin = buildin_default;
end
if (nargin >= 4) & ~isempty(Length_Main)
    all = Length_Main;
else
    all = all_default;
end
if (nargin >= 5) & ~isempty(NumberOfHypotheses)
    MK = NumberOfHypotheses;
else
    MK = MK_default;
end
if (nargin >= 6) & ~isempty(DisplayMode)
    DPMode = DisplayMode;
else
    DPMode = DPMode_default;
end

% Setting hyperparameter alpha and beta
display('Generate system hyperparameters!');
LL = 3000;
PD = 5;
for i=1:LL
    sign = 1;
    while sign
        k = ceil(rand(1,2)*N);
        if k(2)-k(1) >= PD
            index = [k(1):k(2)];
            xx = Y(index);
            mx(i) = mean(xx);
            sign = 0;
        end
    end
end
beta = mean(mx)/var(mx);
alpha = mean(mx)*beta;

% Initializing Output Variables
PM = zeros(MK,1);                               % Prob(Hypothesis)
WM = zeros(1,all);                              % Accepted hypothesis index in each iteration of the main period
SK = zeros(MK-1,all,MK);                        % Samples of change-point index 
SL = zeros(MK,all,MK);                          % Samples of Rate

% Setting System Variable
H = zeros(1,N);                                 % cumulative counts of the series Y; H(i) = sum( Y(j) | j < i)
H(1) = Y(1);
for i=2:N
    H(i) = H(i-1)+Y(i);
end
H = [0,H];

% Initializing Simulation Variables
tmpK = zeros(MK,MK+1);                          % temperory samples of change-point index within an iteration
tmpLBD = zeros(MK,MK);                          % temperory samples of rates within an iteration
tmpK(1,1) = 1;
tmpK(1,2) = N+1;
tmpLBD(1,1) = gamrnd(H(N+1)+alpha,1/(N+beta));
kOld = 1;

% MCMC Simulation Loop
display('The burn-In period has started!');
for i=1:buildin+all
    if (i == buildin + 1)
        display('The main period has started!');
    end
    % Inter-Model move probability
    if kOld == 1
        kNew = 2;
        rGo = 1;
        rBack = 1/2;
        if kNew == MK
            rBack = 1;
        end
    end
    if kOld == MK
        kNew = MK-1;
        rGo = 1;
        rBack = 1/2;
        if kNew == 1
            rBack = 1;
        end
    end
    if kOld > 1 & kOld < MK
        if (rand(1)>0.5) 
            kNew = kOld+1;
        else
            kNew = kOld-1;
        end
        rGo = 1/2;
        rBack = 1/2;
        if kNew == 1 | kNew == MK
            rBack = 1;
        end
    end
    if kNew > kOld 
        % In case of a move to increase dimension by 1
        sign = 1;
        while (sign == 1)
            tmp = ceil(rand(1)*N);
            k = sum(tmp>=tmpK(kOld,1:kOld+1));
            k0 = tmpK(kOld,k);
            k1 = tmpK(kOld,k+1);
            lbd = gamrnd(H(k1)-H(k0)+alpha,1/(k1-k0+beta));
            index = [k0+1:k1-1]-k0; 
            LN = length(index);
            if LN >= 1
                sign = 0;
            end
        end
        M = zeros(1,LN);
        % get k | lbd
        for j=1:LN
            lbd1 = (H(k0+j)-H(k0)+alpha)/(j+beta);
            lbd2 = (H(k1)-H(k0+j)+alpha)/(LN+1-j+beta);
            M(j) = index(j)*(lbd2-lbd1)+log(lbd1/lbd2)*(H(k0+j)-H(k0));
        end
        M = exp(M-max(M));
        M = M/sum(M);
        RM = DisSample(index,M,1);
        PO = log(M(RM));
        PR = log(1/(N-kOld));
        RM = RM+k0;
        % get lbd | k
        lbd1 = gamrnd(H(RM)-H(k0)+alpha,1/(RM-k0+beta));
        lbd2 = gamrnd(H(k1)-H(RM)+alpha,1/(k1-RM+beta));
        PR = PR + log(gampdf(lbd1,alpha,1/beta)+realmin);
        PR = PR + log(gampdf(lbd2,alpha,1/beta)+realmin);
        PR = PR - log(gampdf(lbd,alpha,1/beta)+realmin);
        PO = PO + log(gampdf(lbd1,alpha+H(RM)-H(k0),1/(RM-k0+beta))+realmin);
        PO = PO + log(gampdf(lbd2,alpha+H(k1)-H(RM),1/(k1-RM+beta))+realmin);
        PO = PO - log(gampdf(lbd,alpha+H(k1)-H(k0),1/(k1-k0+beta))+realmin);
        LK = (H(RM)-H(k0))*log(lbd1) - (RM-k0)*lbd1;
        LK = LK + (H(k1)-H(RM))*log(lbd2) - (k1-RM)*lbd2;
        LK = LK - (H(k1)-H(k0))*log(lbd) + (k1-k0)*lbd;;
    else 
        % In case of a move to decrease dimension by 1
        tmp = ceil(rand(1)*N);
        tmp = abs(tmp-tmpK(kOld,2:kOld));
        [A,B] = min(tmp);
        k = B+1;
        k0 = tmpK(kOld,k-1);
        k1 = tmpK(kOld,k+1);
        index = [k0+1:k1-1]-k0; 
        LN = length(index);
        M = zeros(1,LN);
        % get k | lbd
        for j=1:LN
            lbd1 = (H(k0+j)-H(k0)+alpha)/(j+beta);
            lbd2 = (H(k1)-H(k0+j)+alpha)/(LN+1-j+beta);
            M(j) = index(j)*(lbd2-lbd1)+log(lbd1/lbd2)*(H(k0+j)-H(k0));
        end
        M = exp(M-max(M));
        M = M/sum(M);
        RM = DisSample(index,M,1);
        PO = -log(M(RM));
        PR = -log(1/(N-kNew));
        RM = RM+k0;
        % get lbd | k
        lbd1 = gamrnd(H(RM)-H(k0)+alpha,1/(RM-k0+beta));
        lbd2 = gamrnd(H(k1)-H(RM)+alpha,1/(k1-RM+beta));
        lbd = gamrnd(H(k1)-H(k0)+alpha,1/(k1-k0+beta));
        PR = PR - log(gampdf(lbd1,alpha,1/beta)+realmin);
        PR = PR - log(gampdf(lbd2,alpha,1/beta)+realmin);
        PR = PR + log(gampdf(lbd,alpha,1/beta)+realmin);
        PO = PO - log(gampdf(lbd1,alpha+H(RM)-H(k0),1/(RM-k0+beta))+realmin);
        PO = PO - log(gampdf(lbd2,alpha+H(k1)-H(RM),1/(k1-RM+beta))+realmin);
        PO = PO + log(gampdf(lbd,alpha+H(k1)-H(k0),1/(k1-k0+beta))+realmin);
        LK = -(H(RM)-H(k0))*log(lbd1) + (RM-k0)*lbd1;
        LK = LK - (H(k1)-H(RM))*log(lbd2) + (k1-RM)*lbd2;
        LK = LK + (H(k1)-H(k0))*log(lbd) - (k1-k0)*lbd;
    end
    r = PR+LK-PO+log(rBack)-log(rGo);
    if r >= 0
        if kNew > kOld 
            tmpK(kNew,1:k) = tmpK(kOld,1:k);
            tmpK(kNew,k+1) = RM;
            tmpK(kNew,k+2:kNew+1) = tmpK(kOld,k+1:kOld+1);
            tmpLBD(kNew,1:k-1) = tmpLBD(kOld,1:k-1);
            tmpLBD(kNew,k) = lbd1;
            tmpLBD(kNew,k+1) = lbd2;
            tmpLBD(kNew,k+2:kNew) = tmpLBD(kOld,k+1:kOld);
        else
            tmpK(kNew,1:k-1) = tmpK(kOld,1:k-1);
            tmpK(kNew,k:kNew+1) = tmpK(kOld,k+1:kOld+1);
            tmpLBD(kNew,1:k-2) = tmpLBD(kOld,1:k-2);
            tmpLBD(kNew,k-1) = lbd;
            tmpLBD(kNew,k:kNew) = tmpLBD(kOld,k+1:kOld);
        end
        kOld = kNew;  
    else

        r = exp(r);
        if r>rand(1)
            if kNew > kOld 
                tmpK(kNew,1:k) = tmpK(kOld,1:k);
                tmpK(kNew,k+1) = RM;
                tmpK(kNew,k+2:kNew+1) = tmpK(kOld,k+1:kOld+1);
                tmpLBD(kNew,1:k-1) = tmpLBD(kOld,1:k-1);
                tmpLBD(kNew,k) = lbd1;
                tmpLBD(kNew,k+1) = lbd2;
                tmpLBD(kNew,k+2:kNew) = tmpLBD(kOld,k+1:kOld);
            else
                tmpK(kNew,1:k-1) = tmpK(kOld,1:k-1);
                tmpK(kNew,k:kNew+1) = tmpK(kOld,k+1:kOld+1);
                tmpLBD(kNew,1:k-2) = tmpLBD(kOld,1:k-2);
                tmpLBD(kNew,k-1) = lbd;
                tmpLBD(kNew,k:kNew) = tmpLBD(kOld,k+1:kOld);
            end
            kOld = kNew;
        end
    end
    % Discarding burn-in period and record the output within the main period
    if (i > buildin)
        WM(i-buildin) = kOld;
        k = kOld;
        SK(1:k-1,i-buildin,k) = tmpK(k,2:k)';
        SL(1:k,i-buildin,k) = tmpLBD(k,1:k)';
    end

    if (mod(i,1000) == 0 & DPMode == 1)
        Number_Of_Iterations = i;
        display(Number_Of_Iterations);
    end
end

display('The main period is end!');
for i=1:MK
    PM(i) = mean(WM==i);
end

base = 1 / N;
if (DPMode == 1)
    figure;
    bar([0:MK-1],PM,'k');
    axis([-0.5, MK-0.5, 0, ceil(max(PM)/0.1)*0.1+0.1]);
    xlabel('Number of Change-points','fontweight','bold');
    title('Posterior probability of each candidate hypothesis','fontweight','bold');
    [AA,BB] = max(PM);
    if (BB > 1)
        zz = SK(1:BB-1,:,BB);
        zz = zz(:, zz(1,:)~=0);
        panelSetting = [1, 1];
        if (BB == 3)
            panelSetting = [2, 1];
        end
        if (BB == 4 | BB == 5)
            panelSetting = [2, 2];
        end
        if (BB == 6 | BB == 7)
            panelSetting = [2, 3];
        end
        if (BB >= 8)
            p = floor(sqrt(BB - 1))+ 1;
            panelSetting = [p, p];
        end
        figure
        for p = 1 : BB-1
            [nrow, ncol] = size(zz);
            [CC,DD] = hist(zz(p,:),[1:N]);
            subplot(panelSetting(1), panelSetting(2), p);
            bar(timeIndex, CC/ncol, 'k');
            axis([min(timeIndex)-0.5, max(timeIndex)+0.5, 0, ceil(max(CC/ncol)/base)*base+base]);
            xlabel('Time Index','fontweight','bold');
            ylabel('Posterior PMF','fontweight','bold');
        end
    end
end

ProbHypothesis = PM;
SampleOfHypothesisIndex = WM;
SampleOfChangePoint = SK;
SampleOfRate = SL;
return 

% Complementary function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = DisSample(k,PMF,l)

% sample from PMF(k) with length l, where k is the defined region (a vector) and PMF is the realtive PMF 
y = zeros(l,1);
K = length(PMF);
C = zeros(1,K);
C(1) = PMF(1);
temp = rand(l,1);
for i=2:K
    C(i) = C(i-1)+PMF(i);
end
for i=1:l
    l = sum(temp(i) > C)+1;
    if l > K
        l = K;
    end
    y(i) = k(l);
end
return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%