%% Bootstrapping -- classic example from Efron and Tibshirani book ("An Introduction to the Bootstrap")

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
nSamples = 1e5; %number of times to resample (you can play around with this value and see how stable your estimates are)
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
lb = sortedRatio(lowerBoundIndex); %Find actual lower bound value
ub = sortedRatio(upperBoundIndex); %Find actual higher bound value

%Plot histogram
nBins = 50;
figure
histogram(ratio,nBins);
xlabel('Ratio')
%Add lines representing mean and CI of this distribution
meanRatio=line([mean(ratio) mean(ratio)], [min(ylim) max(ylim)],'color','--k','linewidth',2)
CI_lb=line([lb lb], [min(ylim) max(ylim)],'color','k','linewidth',2)
CI_ub=line([ub ub], [min(ylim) max(ylim)],'color','k','linewidth',2)
title('Bootstrapped distribution of HA ratios (Aspirin/Placebo)')

%What are other times you would use bootstrapping? 
%It's great for estimating CIs because it makes no assumptions about the
%distribution of the data (ie, you can get reliable bootstrapped CIs for 
%all types of distributions--skewed, uniform, whatever). Why? think Central
%Limit Theorem. You can apply this idea to statistical tests as well,
%especially when your data are non-normally distributed. We'll next go over
%non-parametric permutation tests, which won't seem as abstract now that
%you understand the basic idea of resampling.