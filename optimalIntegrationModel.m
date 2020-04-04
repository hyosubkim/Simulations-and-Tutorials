%%% Here we simulate predictions of how multisensory integration may 
%%% operate in healthy subjects. This Bayesian cue combination model is
%%% taken from Alais & Barr 2004, although several other studies use a
%%% similar model (e.g., Ernst & Banks 2002, van Beers et al 1996, etc.). 

%Clear your workspace
clear all; close all; clc

%Defining the range of probes
res = 0.1; %pretending we're testing at every 1/10 of a degree (not actually practical, though)
stim = -20:res:20;

%We're simulating a single trial here 
midline = 0; %center of workspace
delta = 5; %how many degrees rightward is visual stimulus presented on conflict presentation
muV = midline+delta; %We assume likelihoods are centered on actual stim location 
muA = midline-delta; 

%These represent the experimentally obtained variances, based on unimodal
%psychometric testing. 
sigmaV = 36;
sigmaA = 6;

%corresponding likelihood fxns -- note that we don't need to explicitly add 
%normal constant term (1/sqrt(2*pi*sigma^2)) bc fxns need to be normalized 
%anyway
pV = exp(-(stim-muV).^2/(2*sigmaV^2));
pV = pV/sum(pV); %normalize probability distribution so that integral is 1
pA = exp(-(stim-muA).^2/(2*sigmaA^2));
pA = pA/sum(pA); %normalize probability distribution so that integral is 1

%calculate combined likelihood function, p(xA,xV|s) -- Note that this is
%equivalent to posterior bc there is no prior probability (ie, all stim
%locations are equi-probable; more explicitly: p(s|xA,xV)=p(xA|s)*p(xV|s)
%and since we have conditional independence p(s|xA,xV) = p(xA,xV|s)
posterior = pA.*pV;  
posterior = posterior/sum(posterior); 

%sanity check: 
wA = sigmaV^2/(sigmaA^2 + sigmaV^2);
wV = sigmaA^2/(sigmaA^2 + sigmaV^2);
s_hat = wA*muA + wV*muV

%determine maximum a posterior estimate (i.e., mean of posterior)
map = sum(stim.*posterior)

%plot likelihoods and posteriors
figure; hold on
plot(stim,pV,'b','linewidth',3)
plot(stim,pA,'r','linewidth',3)
plot(stim,posterior,'k','linewidth',3)
plot([map map],[min(ylim) max(ylim)], '--k')
xlabel(['Probe displacement (',char(176), ')'])
ylabel('P')
legend({'likelihood V', 'likelihood P', 'posterior', 'MAP'})
legend('boxoff')
title('Optimal integration')


