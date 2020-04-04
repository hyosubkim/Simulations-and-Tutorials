%%%%% Relevance estimation model of Wei and Kording (2009)%%%%%

%This script was adapted from code provided by Konrad Kording at the
%Computational Sensory-Motor Neuroscience Summer School (2017). All the
%amazing course materials are freely available here:
%http://compneurosci.com/wiki/index.php/CoSMo_2017

clear all; close all; clc

%all the values that the subject could have perceived
x=-20:0.1:20; %we always need to integrate over unobserved signal in brain
perturbations=[-8    -4    -2    -1     0     1     2     4     8  ];


%Proprioceptive likelihood
muP=0;
sigmaP=10; %does not change anything as only ratios matter
pP=exp(-(x-muP).^2/(2*sigmaP^2));
pP=pP/sum(pP);


%Loop over stimulus parameters (visual perturbation)
for i=1:length(perturbations)
    
    figure; hold on
    p1=plot(pP,'r','linewidth',3);
    xlim([0 length(x)+1])
    set(gca,'xtick',1:40:length(x),'xticklabel',-20:4:20)
    xlabel('Possible hand locations (cm)')
    ylabel('P')
    
    %For visual likelihood
    sigmaV=5;
    muV=perturbations(i); %visual likelihood is centered on perturbed visual feedback location
    pV=exp(-(x-muV).^2/(2*sigmaV^2));  %this is Gaussian apart from a constant
    pV=pV/sum(pV);  % probabilities must add up to 1
    p2=plot(pV,'g','linewidth',3); hold on;
    
    %combined variance
    sigmaCombined = sqrt((sigmaV^2*sigmaP^2)/(sigmaV^2+sigmaP^2));
  
    %p(relevant) aka p(causal)
    C=1; %turns out this has effectively the same effect as third para
    pRelevant(i,1)=exp(-(muV).^2/(2*sigmaCombined^2))./(exp(-(muV).^2/(2*sigmaCombined^2))+C);
    
    %combine relevant and irrelevant cues
    xHat=pRelevant(i)*(pV.*pP)/sum(pV.*pP)+(1-pRelevant(i))*pP;
    p3=plot(xHat,'k','linewidth',3); % all the estimates that the subject could have seen
    meanPosterior(i)=sum(xHat.*x);
    scaling=-.5; %characterizes the magnitude of the influence of a visual disturbance on future trials
    predictedCorrection(i)=scaling*meanPosterior(i);
    legend({'p(xProp)','p(xVis)','xHat'})
    legend('boxoff')
end


figure; hold on
plot(perturbations, predictedCorrection,'b','linewidth',3);
xlabel('Perturbations (cm)')
ylabel('Corrections (cm)')

figure; hold on
plot(perturbations,pRelevant./max(pRelevant),'k','linewidth',3)
xlim([-10 10])
xlabel('Perturbations')
ylabel('Normalized p(relevant)')
