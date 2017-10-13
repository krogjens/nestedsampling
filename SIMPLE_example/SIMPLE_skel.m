function SIMPLE_skel()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This skeleton routine carries out Bayesian inference 
% for the data specified in
% misc.data_id together with data path. A number of parameters
% may be specified by the user:
%
%  - Dmin, Dmax, mumin, mumax
%    specify the minimum and maximum values of the two parameters   
%    which model the dynamics, and must be chosen to ensure that 
%    the most likely volume of parameter space is within these ranges.
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
%  - the "models" structure specify the functions used for likelihood calculation,
%    sample generation, which are specific 
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

% There are two models:
% Model 1 is Brownian motion
% Model 2 is Brownian motion with a drift

%Specify a Jeffreys prior on the diffusion constant
Dmin=10^(-1);
Dmax=10^2;   
models(1).invprior=@(u) Dmin*exp(u*log(Dmax/Dmin));

%Additionally specify a uniform prior on the drift for Model 2
mumin=-10;
mumax=10;
models(2).invprior=@(u) [Dmin*exp(u(1)*log(Dmax/Dmin)); (mumax-mumin)*u(2)+mumin];

%Specify options
options.nwalkers=200;   % Number of walkers to be generated
options.stoprat=10^(-3);% Ratio to stop sampling
options.nsteps=30;      % Attempted number of steps in parameter space
models(1).options=options;
models(2).options=options;

%Specify the u-generators
models(1).genu=@() rand(1,1);
models(2).genu=@() rand(1,2);

%Specify the logl and logl_n
log_normal=@(obs,mu,var) -log(sqrt(2*pi*var))*(length(obs)-1)-sum((obs(2:end)-obs(1:(end-1))-mu).^2)/(2*var);
models(1).logl=@(obs,theta) log_normal(obs,0,2*theta(1));
models(1).logl_n=@(obs,theta,n) log_normal(obs,0,2*theta(1)*n);
models(2).logl=@(obs,theta) log_normal(obs,theta(2),2*theta(1));
models(2).logl_n=@(obs,theta,n) log_normal(obs,theta(2)*n,2*theta(1)*n);

%Specify the index for the labels
models(1).labels=[1];
models(2).labels=[1 2];

%Percentiles to be calculated
misc.percentiles_at=[0.02 0.16 0.5 0.84 0.98];

%Labels for the parameters and percentile line in the output text-file
misc.labels=...
['Diffusion constant: ';...
 'Drift:              '];
misc.titles=...
{'Percentile at:      '};

%Specify output filename beginnings
misc.data_id = 'simple';
data=cumsum(sqrt(2*1)*randn(1000,1));

%Tell ns_print to write a summary-file
misc.nssummary=['_results.txt'];

%Run the nested sampling algorithm for all models and compile results
[results] = ns_processdataset(data,models,misc);

path=[misc.data_id,'_output'];

%Save results
save(path,'results')
save(path,'data','-append')
