%% Bootstrapping -- classic example from Efron and Tibshirani book ("An Introduction to the Bootstrap")

clear all; close all; clc

% This example has been adapted from an excellent series of lectures that 
% goes in depth with the theory and application of resampling methods and 
% much much more: http://www.cns.nyu.edu/~eero/math-tools/. 

% Q: Does taking aspirin reduce your risk of heart attack? 
% These are data from a famous study done back in the 80s.
nA = 11037; %total number of participants taking aspirin
nAH = 104; %number of aspirin takers who ended up having a heart attack
nP = 11034; %total number of participants on a placebo
nPH = 189; %number of placebo takers who had a heart attack
empiricalRatio = (nAH/nA)/(nPH/nP); %experimentally observed ratio of HAs in Aspirin vs Placebo groups (i.e., Experimental Results)
nSamples = 1e4; %number of times to resample (you can play around with this value and see how stable your estimates are)
ratio = zeros(nSamples,1); %initialize and preallocate (makes code run faster)

A1 = ones(nAH,1); %we code those who had a heart attack with a 1
A2 = zeros(nA-nAH,1);%aspirin takers who did not have a heart attack
P1 = ones(nPH,1); %same for placebo
P2 = zeros(nP-nPH,1); %same 

%Stack your vectors of aspirin and placebo groups
A = [A1;A2]; 
P = [P1;P2]; 


%Here is where the magic takes place. We're going to simulate the same
%experiment a total of nSamples times by randomly resampling with
%replacement, thereby creating the distribution of sample means (Question:
%Do we expect a normal distribution? Why or why not?)
for ii = 1:nSamples 
    indicesA = randi(nA,[nA 1]); %draw random integers with max value equal to number of observations in sample
    indicesP = randi(nP,[nP 1]); %same idea
    ratio(ii) = sum(A(indicesA))/sum(P(indicesP)); %ratio of HAs in simulated sample
end

%We now have a bootstrapped distribution. Let's calculate CIs.
%First step: Sort values
lbPercentile = 2.5; ubPercentile = 97.5; %you could change these if you want something other than 95% CIs
sortedRatio = sort(ratio); 
lbIdx = round(nSamples/100*lbPercentile); %Find index of the lower bound
ubIdx = round(nSamples/100*ubPercentile); %Find index of upper bound 
lb = sortedRatio(lbIdx); %Find actual lower bound value
ub = sortedRatio(ubIdx); %Find actual higher bound value

%Plot histogram
nBins = 50;
figure
histogram(ratio,nBins);
xlabel('Ratio')
%Add lines representing mean and CI of this distribution
meanRatio=line([mean(ratio) mean(ratio)], [min(ylim) max(ylim)],'color','k','linestyle','--','linewidth',2)
CI_lb=line([lb lb], [min(ylim) max(ylim)],'color','k','linewidth',2)
CI_ub=line([ub ub], [min(ylim) max(ylim)],'color','k','linewidth',2)
title('Bootstrapped distribution of HA ratios (Aspirin/Placebo)','fontsize',14)

%%%What are other times you would use bootstrapping?%%% 
%It's great for estimating CIs because it makes no assumptions about the
%distribution of the data (ie, you can get reliable bootstrapped CIs for 
%all types of distributions--skewed, uniform, whatever). Why? think Central
%Limit Theorem. You can apply this idea to statistical tests as well,
%especially when your data are non-normally distributed. We'll next go over
%non-parametric permutation tests, which won't seem as abstract now that
%you understand the basic idea of resampling.

%%%Where have we seen bootstrapping already?%%%
%Taylor & Ivry 2011 (aiming model); Smith et al. 2006 (dual-rate model);
%Charalambous et al. 2018 (fitting behavioral data from chronic stroke
%participants)

%%%Sample applications of bootstrapping%%%
%1) You want to calculate confidence intervals on a sample statistic coming
%from some distribution you are not sure about
%2) You performed a correlational anlaysis, and you want to get some sense
%of the uncertainty of your correlation coefficient, or regression weights
%3) You want to figure out CIs on parameter estimates from a model--e.g.,
%you fit a state-space model, where A and B are constrained to be between 0
%and 1
%4) Just about anything else you can think of...


%%%If bored, see also the Jupyter notebook version of this tutorial%%%

%% Permutation tests

clear all; close all; clc

%What if you have non-parametric data and want to compare the difference
%between two groups? Or if you want to compare three or more groups? How
%about if you have run a regression analysis and want to whether the linear
%relationship between variables is reliable? This is what permutation tests
%are for--and they're not restricted to non-parametric data...Permutation 
%tests are distribution free. We use the data itself to approximate the 
%distribution of sample whatevers. 

%Here we'll look at BDNF data from Charalambous et al. 2017. 
%Read in data from Excel spreadsheet
[status,sheets] = xlsfinfo('PrimReSt_SNP_X_Learning.xlsx');
T=readtable('PrimReSt_SNP_X_Learning.xlsx','Sheet',1);

%first test on Controls
group = categorical(T.group);
genotype = categorical(T.genotype);

CV = T.changeOfRate(group=='CON' & genotype=='VAL')'; %change of adaptation rate for ctls with Val
CM = T.changeOfRate(group=='CON' & genotype=='MET')'; %change of adaptation rate for ctls with Met

TV = T.changeOfRate((group=='TMW' | group=='TBE') & genotype=='VAL')'; %change of adaptation rate for treatment grp VALs
TM = T.changeOfRate((group=='TMW' | group=='TBE') & genotype=='MET')'; %change of adaptation rate for treatment grp METs


numMCS = 1e5; %Number of Monte Carlo simulations (try increasing to at least 1e6 and see how stable/unstable p-vals are.
              %The name Monte Carlo comes from the famous gambling resort area--not to be confused with Markov Chain Monte
              %Carlo sampling methods.)
nCV = length(CV); %Number of ctl Vals
nCM = length(CM); %Number of ctl Mets

nTV = length(TV); %Number of tx VALs
nTM = length(TM); %Number of tx METs

%Step 1 in this process: Come up with a reasonable test statistic to
%compare (difference in groups seems reasonable)
testStatistic1 = mean(CV)-mean(CM)
testStatistic2 = mean(TV)-mean(TM)

%A: How likely is a value of xx (the value we got) due to chance alone? We
%permute the values and draw groups arbitrarily (randomly) and see what
%values of the test statistic we get. Doing this a large number of times
%allows us to create a null distribution of the test statistic. 
%Logic: If there is no effect of the genotype, it should not matter what 
%group the adaptationRates came from.

dTS1 = zeros(numMCS,1); %Initialize and preallocate the distribtuion of the test staistic
dTS2 = zeros(numMCS,1);
allData1 = [CV CM]; %Concatenate all the ctl data
allData2 = [TV TM]; %Concatenate all the tx data

%Monte Carlo method: Random permutation/reshuffling of the group labels;
%like bootstrapping but here we are sampling WITHOUT replacement
for ii=1:numMCS %However many times you want to do this
    indices1 = randperm(nCV + nCM); 
    indices2 = randperm(nTV + nTM);
 
    tempDATA1 = allData1(indices1); %Create a newly permuted dataset. We now parse it into groups.
    group1 = tempDATA1(1:nCV); %Capture the first n elements
    tempDATA1(:,1:nCV) = []; %Delete those from the permuted dataset
    group2 = tempDATA1; %Group 2 is simply what is left 
    
    tempDATA2 = allData2(indices2);
    txgroup1 = tempDATA2(1:nTV);
    tempDATA2(:,1:nTV) = [];
    txgroup2 = tempDATA2;
    
    dTS1(ii) = mean(group1)-mean(group2); %capture the simulated test statistic for each run
    dTS2(ii) = mean(txgroup1)-mean(txgroup2);
end

figure
histogram(dTS1,40)

%Finding the "exact p-value"
%exactP1 =sum((dTS1>=testStatistic1))/length(dTS1); %One-tailed
%Out of curiosity: Two-tailed case
exactP1 = sum(testStatistic1<=abs(dTS1))/length(dTS1);

%Add some labels
line([testStatistic1 testStatistic1],[min(ylim) max(ylim)],'color','k','linewidth',3)
line([-testStatistic1 -testStatistic1],[min(ylim) max(ylim)],'color','k','linewidth',3)
xlabel('Distribution of simulated test statistics in absence of effect')
ylabel('Frequency')
title(['Distribution of test statistics (Controls). Exact p = ', num2str(exactP1,'%1.5f')])


figure
histogram(dTS2,40)

%Finding the exact p value
%exactP2 =sum((dTS2>=testStatistic2))/length(dTS2); %One-tailed
%Out of curiosity: Two-tailed case
exactP2 = sum(testStatistic2<=abs(dTS2))/length(dTS2);

%Add some labels
line([testStatistic2 testStatistic2],[min(ylim) max(ylim)],'color','k','linewidth',3)
line([-testStatistic2 -testStatistic2],[min(ylim) max(ylim)],'color','k','linewidth',3)
xlabel('Distribution of simulated test statistics in absence of effect')
ylabel('Frequency')
title(['Distribution of test statistics (Tx Groups). Exact p = ', num2str(exactP2,'%1.5f')])


