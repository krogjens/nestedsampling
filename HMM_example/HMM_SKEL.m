function HMM_SKEL()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hidden Markov Model for different diffusive states (brownian motion)
%Parameters of a model are 
%Diffusion constant (D) for each state
%Drift constant (mu) for each non-zero drift state
%Transition probability matrix (T)
%
%Different models are distinquished by total number of states, and number
%of non-zero drift states 
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The ordering of the parameters in the list will be 
% 1. n diffusion parameters
% 2. j drift parameters
% 3. state transition parameters, which are indices of a transition matrix 
% BUt only the off diagonals are active variables
%and finally #states  starting probsbility variables (i think necesarry
%because uniform prior)
% example of 3-state, 2-drift hmm
%thetas = [sigma1,sigma2,sigma3,mu1x,mu1y,mu2x,mu2y,T_11,T_21,T_31,T_12,T_22,T_32,T_13,T_23,T_33,P_1,P_2,P_3]
%By passing along the number of states, and drift containing states, we can
%manage all the parameters as a list.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This skeleton routine carries out Bayesian inference and
% information content model checks for the data specified in
% misc.data_id together with data path. A number of parameters
% may be specified by the user:
%
%  - Dmin, Dmax, amin, amax, fmin, fmax, noisemin, and noisemax
%    specify the minimum and maximum values of the four parameters   
%    which model the dynamics, and must be chosen to ensure that 
%    the most likely volume of parameter space is within these ranges.
%
%  - Tmin is the minimum length of a rescaled dataset to be used for
%    for the information content check.
%
%  - options.nwalkers is the number of initial samples used in the
%    nested sampling algorithm to sample parameter space.
%
%  - options.stoprat is the ratio between the last added evidence and the
%    total samples evidence at which nested sampling will terminate.
%
%  - options.nsteps is the number steps attempted with the MCMC
%    technique to generate a uniformly distributed sample from another sample.
%
%  - nlist is the list of rescalings used to carry out the information
%    content check.
%
%  - trackmax is the number of tracks replicated for each scaling in the
%    in the information content check.
%
%  - the "models" structure specify the functions used for likelihood calculation,
%    sample generation, data rescaling, and data replication, which are specific 
%    to the system of interest. To use the nested sampling framework on another system, 
%    similar functions must exist and be specified in this structure.
%
% The routine outputs a .mat file with information on the user set parameters
% as well as evidence estimations and inferred parameter means and the samples used.
% In addition a .txt file is written, holding the conclusions of the analysis.
% The routine writing this file takes the following inputs:
%  - misc.percentiles are the percentiles used for the characterization of the
%    posterier 
%  - misc.labels are the labels assigned to each inferred parameter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%Specify location of data
misc.data_id = 'data/hmmdata';

data_path = [misc.data_id,'.csv'];

%load data 

%format is a series of x and y displacements
data = importdata(data_path);
data=data.data;
%data = data.obs; % Uncomment this if running on noisetrack.mat or
                  % largenoisetrack.mat
realA=[0.95,0.05;0.03,0.97];
realp=[0.5;0.5];
realB=[0.8;0.4];
realu=[0,0;0,0];
nsteps=500;
dt=1/10; % 1/sampling rate in hz
data=zeros(nsteps,2);            
[data(:,1),data(:,2),~]= HMMgenerator(realp,realA,realB,realu,nsteps);
%Specify prior ranges
%Diffusion constant
sigmamin=10^(-4);
sigmamax=10^2;   

%Drift constant
mu_min = -2;
mu_max = 2;


%Convert prior ranges to dimensions of pixels and frames for convenience 
% ranges - a 2x2 array of minimum and maximum values for the 2 parameters:
%  - standard deviation of gaussian, and the mean mu
ranges=[sigmamin sigmamax ; mu_min mu_max];

Tmin = 10; % Required minimum length of a scaled data set
nmax = floor((length(data)-1) / Tmin);
if ~(nmax > 1)
   fprintf('Data set too short to be scaled to length >= %i',Tmin);
elseif nmax > 10
   nmax = 10;
end


%Specify options
options.nwalkers=300;   % Number of walkers to be generated
options.stoprat=10^(-6);% Ratio to stop sampling
options.nsteps=30;      % Attempted number of steps in parameter space
options.nlist = 1:1:nmax;    % Maximum time scaling factor for data
options.trackmax = 100; % Number of replicated tracks to compare information for
options.Nparfor=6;
final=@(x) x(end); % extract last element in a vector

%Max number of states to fit
Jtot = 2; 

%Specify the models
models=repmat(struct('states',0,'drifts',0),1,(Jtot^2+3*Jtot)/2);
%create selection which tells the number of states, and how many with drift
for M=1:Jtot
    for d=0:M
        models(0.5*M^2+0.5*M+d).states=M;
        models(0.5*M^2+0.5*M+d).drifts=d;
    end
end

for i=1:(Jtot^2+3*Jtot)/2
  models(i).options=options;
  
  %number of u values that will be generated. (diagonal value in transition
  %matrix will be generated also, but thrown out. This is for the sake of
  %ease of implementation and speed
  %nstates diffusion parameters
  %ndrifts drift parameters
  %nstates^2-kroneckerdelta(1,nstates) transition rates
  %nstates-kroneckerdelta(1,nstates) starting probabilities
  models(i).genu=@() HMM_generate_u(models(i).states+models(i).drifts*2+models(i).states^2-kroneckerdelta(1,models(i).states)+models(i).states-kroneckerdelta(1,models(i).states));
  %logl function, is the forward algorithm for multiple states or regular gaussian for 1 state. C code is preferable as it is
  %terribly slow as a .m 
  models(i).logl=@(obs,theta) HMM_logl(obs,theta, models(i).states, models(i).drifts, length(data),dt);
  
  %inverse prior for HMM, we use
  models(i).invprior=@(u) HMM_invprior(u,ranges,models(i).states,models(i).drifts);
  %label the models as #states #drift
  models(i).labels=[models(i).states:models(i).drifts];
end

%Percentiles to be calculated
misc.percentiles_at=[0.02 0.16 0.5 0.84 0.98];

%Labels for the parameters
misc.labels=...
['Step deviation: ';...
 'Exponent:       ';...
 'Force:          ';...
 'Measurem. err.: '];

misc.titles=...
{'Percentile at:  '};

%Tell ns_print to write a summary-file
misc.nssummary=['_results.txt'];

%Run the nested sampling algorithm for all models and compile results
[results] = ns_processdataset(data,models,misc);

path=[misc.data_id,'_output'];

%Save results
save(path,'results')
save(path,'data','-append')
save(path,'options','-append')
save(path,'ranges','-append')
