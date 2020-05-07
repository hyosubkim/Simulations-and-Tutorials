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


clear all; close all; clc

load('adaptationData.mat')

%% Plot group analyzed hand data
zeroline = zeros(size(hand_ang_m,1));
figure(1)
hold on
set(gcf,'units','inches','pos',[5 5 4.75 2.63]);
set(gcf,'PaperPositionMode','auto')
rectangle('position',[0 -10 10 40],'facecolor',gray,'edgecolor',gray)
rectangle('position',[240 -5 10 40],'facecolor',gray,'edgecolor',gray)
for i=1:2
   
    grpdata(i)=shadedErrorBar(1:length(hand_ang_m),nanmean(hand_ang_m(grp==grpid(i),:)),...
        nanstderr(hand_ang_m(grp==grpid(i),:)),{'.','markersize',8,'color',grpclrs(i,:)},1);
    
end
xlim([-2 272])
ylim([-5 30])
plot(xlim, [0 0],'-k',[10 10], ylim,'--k',[20 20], ylim, '--k',...
    [240 240], ylim, '--k',[250 250], ylim, '--k')
title(['Extended 1.75',char(176),' clamp'])
xlabel('Movement cycle (4 reaches)')
ylabel(['Hand angle (',char(176),')'])
% legend([grpdata(1).mainLine,grpdata(2).mainLine],{grplabel{1},grplabel{2}})
% print('smallTgtOnly','-painters','-dpdf')
%print('e2_extTgt','-painters','-dpdf')

%% Cluster permutation test:

% This is a very conservative test to see if you have differences between
% functions in which data at consecutive time points are not independent 
% of each other (eg, adaptation functions). If you have understood
% resampling, in general, and non-parametric permutation tests,
% specifically, you'll see that this is a logical extrapolation of those
% methods. They use this test in Labruna et al. (2019); to my knowledge,
% this is the first instance of using this method in an adaptation study.
% It's more commonly used in studies incorporating EEG and other forms of 
% neural recordings. 

% Warning: This script takes a long time to run even when resampling < 1000
% times because you're working with complete time series data rather
% than discrete values as in previous examples. 

clear all; close all; clc

%load some data from Kim et al. (2019)
load('adaptationData.mat')

%Plot group analyzed hand data to see two groups adapting to a visual clamp
%perturbation
zeroline = zeros(size(hand_ang_m,1));
figure
hold on
set(gcf,'units','inches','pos',[5 5 4.75 2.63]);
set(gcf,'PaperPositionMode','auto')
rectangle('position',[0 -10 10 40],'facecolor',gray,'edgecolor',gray)
rectangle('position',[240 -5 10 40],'facecolor',gray,'edgecolor',gray)
for i=1:2
   
    grpdata(i)=shadedErrorBar(1:length(hand_ang_m),nanmean(hand_ang_m(grp==grpid(i),:)),...
        nanstderr(hand_ang_m(grp==grpid(i),:)),'lineprops',{'.','markersize',8,'color',grpclrs(i,:)},'transparent',1);
    
end
xlim([-2 272])
ylim([-5 30])
plot(xlim, [0 0],'-k',[10 10], ylim,'--k',[20 20], ylim, '--k',...
    [240 240], ylim, '--k',[250 250], ylim, '--k')
title(['Extended 1.75',char(176),' clamp'])
xlabel('Movement cycle (4 reaches)')
ylabel(['Hand angle (',char(176),')'])
% legend([grpdata(1).mainLine,grpdata(2).mainLine],{grplabel{1},grplabel{2}})


%moving on to the cluster permutation
numBslCycles = 20;
numClampCycles = 220;

%We're going to look at the data cycle-by-cycle and run a continuous series
%of t-tests (we could use other measure, like difference scores, but t-stat
%is more conservative)
for i=1:numClampCycles
    [~,p(i),~,stats] = ttest2(hand_ang_m(grp==grpid(1),numBslCycles+i),hand_ang_m(grp==grpid(2),numBslCycles+i));
    t_score(i) = stats.tstat;   
end
%creating a logical vector to identify clusters of reliable group
%differences
sig_pval = p < .05;

%we pad with 0's at beginning and end to get right indices
sig_pval = [0 sig_pval 0];

%find the clusters
edges = diff(sig_pval);
rising = find(edges==1);
falling = find(edges==-1);

%check the length of each cluster
clusterLength = falling - rising;

%make sure cluster has at least 2 adjacent cycles where group differences
%are reliable
wideEnough = clusterLength >= 2;
startPos = rising(wideEnough);
endPos = falling(wideEnough)-1;  
for j=1:length(startPos)
    %sum the t-scores from each significant cluster to get the "mass" of
    %each cluster
    clusterMass(j) = sum(abs(t_score(startPos(j):endPos(j))));
end
[tsum_data,max_cluster_idx] = max(clusterMass);
sig_cycles = [startPos(max_cluster_idx): endPos(max_cluster_idx)];

%plot boundaries of largest cluster
figure(1)
plot([startPos(max_cluster_idx)+numBslCycles, startPos(max_cluster_idx)+numBslCycles],...
    [min(ylim) max(ylim)],'--m')
plot([endPos(max_cluster_idx)+numBslCycles, endPos(max_cluster_idx)+numBslCycles],...
    [min(ylim) max(ylim)],'--m')    


%plot p-values
figure(2); hold on
plot(p,'k','linewidth',1.5)
plot(xlim,[0.05 0.05],'-r')
xlabel('Cycle number')
ylabel('p-value')


%create the null distribution
nReshuffles = 1e2;
for ii=1:nReshuffles
    
    %shuffle group assignmenets on every iteration
    perm_grp = grp(randperm(length(grp)));
    
    %we're doing the exact same thing as before, but now with shuffled data
    %so that we can create the null distribution
    for jj=1:numClampCycles
        [~,p_null(jj),~,stats_null] = ttest2(hand_ang_m(perm_grp==grpid(1),numBslCycles+jj),hand_ang_m(perm_grp==grpid(2),numBslCycles+jj));
        tscore_null(jj) = stats_null.tstat;
    end
    
    %I should probably rename these variables since I used them before, but
    %I'm leaving them the same for now
    sig_pval = p_null < .05;
    sig_pval = [0 sig_pval 0];
    edges = diff(sig_pval);
    rising = find(edges==1);
    falling = find(edges==-1);
    clusterLength = falling - rising;
    wideEnough = clusterLength >= 2;
    startPos = rising(wideEnough);
    endPos = falling(wideEnough)-1;
    
    if startPos
        atLeastOneCluster(ii) = 1;
        for kk=1:length(startPos)
            t_all_sig = sum(abs(tscore_null(startPos(kk):endPos(kk))));
        end
        tsum_null(ii) = max(t_all_sig);
    else
        tsum_null(ii) = 0;
    end
    
end

exact_pval = sum(tsum_null>tsum_data)/nReshuffles

%plot null distribution
figure; hold on
histogram(tsum_null, 20)
plot([tsum_data tsum_data],ylim,'--r')
xlabel('Cluster masses of null distribution')
ylabel('Frequency')
title(['P-val: ',num2str(exact_pval)])




