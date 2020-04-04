%%%%% Relevance estimation model of Wei and Kording (2009)%%%%%

%This script was adapted from code provided by Konrad Kording at the
%Computational Sensory-Motor Neuroscience Summer School (2017). All the
%amazing course materials are freely available here:
%http://compneurosci.com/wiki/index.php/CoSMo_2017

clear all; close all; clc

%all the values that the subject could have perceived
x=-10:0.1:10; %we always need to integrate over unobserved signal in brain
perturbations=[-8    -4    -2    -1     0     1     2     4     8  ];


%Proprioceptive likelihood
muP=0;
sigmaP=1; %does not change anything as only ratios matter
pP=exp(-(x-muP).^2/(2*sigmaP^2));
pP=pP/sum(pP);

figure; hold on
plot(pP,'r');
xlim([0 length(x)+1])
set(gca,'xtick',1:20:length(x),'xticklabel',-10:2:10)
xlabel('Possible hand locations (cm)')
ylabel('P')

%For visual likelihood
sigmaV=4.25;

%Loop over stimulus parameters (visual perturbation)
for i=1:length(perturbations)
    muV=perturbations(i); %visual likelihood is centered on perturbed visual feedback location
    pV=exp(-(x-muV).^2/(2*sigmaV^2));  %this is Gaussian apart from a constant
    pV=pV/sum(pV);  % probabilities must add up to 1
    plot(pV,'g'); hold on;
  
    %p(relevant) aka p(causal)
    C=1; %turns out this has effectively the same effect as third para
    pRelevant(i,1)=exp(-(muV).^2/(2*sigmaV^2))./(exp(-(muV).^2/(2*sigmaV^2))+C);
    
    %combine relevant and irrelevant cues
    posterior=pRelevant(i)*(pV.*pP)/sum(pV.*pP)+(1-pRelevant(i))*pP;
    plot(posterior,'k'); % all the estimates that the subject could have seen
    meanPosterior(i)=sum(posterior.*x);
    scaling=-5.5535; %characterizes the magnitude of the influence of a visual disturbance on future trials
    prediction(i)=scaling*meanPosterior(i);   
end

figure; hold on
plot(perturbations, prediction,'b');
xlabel('Perturbations (cm)')
ylabel('Corrections (cm)')

figure; hold on
plot(perturbations,pRelevant./max(pRelevant),'k','linewidth',3)
xlim([-10 10])
xlabel('Perturbations')
ylabel('Normalized p(relevant)')
