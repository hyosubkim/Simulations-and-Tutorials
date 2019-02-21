% clear your workspace
clear all; close all; clc

%%% We're going to simulate an experimental session and see if we can
%%% recover the parameters for Tp (movement preparation time) using maximum 
%%% likelihood estimation.

N = 300; %num trials

% We're going to create a subject whose mu and sigma for Tp are 130 and 25,
% respectively. We're actually creating a probability distribution object 
% which we'll pass to pdf.
Tp = makedist('normal',130,25); 

% Prep times are uniformly distributed between 0 and 300 ms and we will 
% assume our participant performs perfectly (i.e.,always moves exactly on 
% the 4th tone.)
RT = randi(300,N,1); 

Tp_pd = pdf(Tp,0:300);    % creates a probability density function
Tp_cd = cumsum(Tp_pd); % creates a cumulative density function

% As always, you want to do a sanity check by plotting. You should see a
% Gaussian (bell-shaped) and a sigmoidal (S-shaped) function
figure(1)
subplot(2,1,1); hold on
real=plot(Tp_pd,'linewidth',2)
ylabel('P')
subplot(2,1,2)
plot(Tp_cd,'linewidth',2)
xlabel('RT (ms)')
ylabel('Cumulative Prob')


% Getting to the model...
alpha = 0.95;   % This is the ceiling level for accuracy. It's constraining
                % max accuracy since it's impossible for people to be 100%
                % accurate.
                
% Let's fit our simulated data
for i=1:N
    % We'll talk about this in journal club
    pH(i) = alpha*(Tp_cd(RT(i)+1)) + (1/8)*(1-Tp_cd(RT(i)+1)); % add 1 to idx 
    H(i,1) = pH(i) > rand;
end
% fminshearch will call your max likelihood function (loglik) and find the
% best fit parameters (remind me to explain fminsearch and why we use
% negative log likelihoods)
[params] = fminsearch(@(x) loglik(x(1),x(2),RT,H), [100,100])   

% plot the fit distributions
Tp_fit = makedist('normal',params(1), params(2));
pdf_Tp_fit = pdf(Tp_fit,0:300);
figure(1); hold on
subplot(2,1,1)
model=plot(pdf_Tp_fit,'g','linewidth',2)
legend([real,model],{'Real','Model fit'})

% How do the parameters look? Is the mean close to 130, and SD close to 25?
% Run the for loop again. Did you get different values? Why?
                
%%

% Let's run the experiment a bunch of times
nSims = 1e2; % set the number of simulated experiments you want to run
for k=1:nSims
    for i=1:N
        % We'll talk about this in journal club
        pH(i) = alpha*(Tp_cd(RT(i)+1)) + (1/8)*(1-Tp_cd(RT(i)+1));
        H(i,1) = pH(i) > rand;
    end
    [params(k,:)] = fminsearch(@(x) loglik(x(1),x(2),RT,H), [100,100]);
end

% Now check the mean of params. Is it pretty close to the actual
% distribution you created? If so, why? 
