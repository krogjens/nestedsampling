function SDE_skel()

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
misc.data_id = 'data/ftrack_data';

data_path = [misc.data_id,'.txt'];

%load data 

data = importdata(data_path);
%data = data.obs; % Uncomment this if running on noisetrack.mat or
                  % largenoisetrack.mat

%Specify prior ranges
Dmin=10^(-4);
Dmax=10^2;   

amin = -2;
amax = 2;

fmin = -1;
fmax = 1;

noisemin=0; 
noisemax=100;

%Convert prior ranges to dimensions of pixels and frames for convenience 
% ranges - a 2x5 array of minimum and maximum values for the 4 parameters:
%   - Mobility coefficient, mobility exponent, force, and measurement error
ranges=[Dmin Dmax ; amin amax; fmin fmax; noisemin noisemax];

Tmin = 10; % Required minimum length of a scaled data set
nmax = floor((length(data)-1) / Tmin);
if ~(nmax > 1)
   fprintf('Data set too short to be scaled to length >= %i',Tmin);
elseif nmax > 10
   nmax = 10;
end

%Specify scalars to calculate p-values
checks=[];

endpos.scalar=@(track) track(end);
endpos.misc.labels={'Model check for final position, p-value: '};
checks=[checks endpos];

nsteps=[1 2 4 8 16];
stepcor.scalar=@(track) util_stepcorrelations(track,nsteps);
misc.labels={'Model check for estimated n-step correlations:'...
;' n:      '... % Note: if this line is deleted then ' Input:  ' is printed
};
misc.columns=nsteps;
stepcor.misc=misc;
checks=[checks stepcor];

nsteps=[1 2 4 8 16];
percentiles=[0.05 0.1:0.1:0.9 0.95];
stepdist.scalar=@(track) util_stepdist(track,nsteps,percentiles);
misc.labels={'Model check with p-values for percentiles (C) for n-step (R) distributions:'...
;' R\\C'... % Note: if this line is deleted then ' R  \  C' is printed
};
misc.columns=percentiles;
misc.rows=nsteps;
stepdist.misc=misc;
checks=[checks stepdist];

%Specify options
options.nwalkers=100;   % Number of walkers to be generated
options.stoprat=10^(-1);% Ratio to stop sampling
options.nsteps=10;      % Attempted number of steps in parameter space
options.nlist = 1:1:nmax;    % Maximum time scaling factor for data
options.trackmax = 100; % Number of replicated tracks to compare information for

final=@(x) x(end); % extract last element in a vector

%Specify the models
MM = [0 0; 1 0; 0 1; 1 1];
n_perm = 2; % Number of permanently active parameters
for i=1:size(MM,1)
  models(i).options=options;
  models(i).genu=@() SDE_generate_u(sum(MM(i,:))+n_perm);
  models(i).logl=@(obs,theta) SDE_logl(obs,SDE_params(theta,MM(i,:)));
  models(i).logl_n=@(obs,theta,n) SDE_logl_m(obs,SDE_params(theta,MM(i,:)),n);
  models(i).scaling =@(obs,n) SDE_scaling(obs,n); % Function that scales data
  models(i).replicate =@(obs,theta,n) SDE_replicate(obs,SDE_params(theta,MM(i,:)),n);
  models(i).checks=checks;
  models(i).invprior=@(u) SDE_invprior(u,ranges,MM(i,:));
  models(i).labels=[1:n_perm];
  for j=1:length(MM(i,:))
    if MM(i,j)==1
      models(i).labels=[models(i).labels j+n_perm];
    end
  end
%  models(i).add{1}=@(theta) theta(1)^2/(2*tau^(2*final(fbm_params(theta,MM(i,:)))));
 % models(i).labels=[models(i).labels 6];
  for j=1:2
    if MM(i,j)==1
%      models(i).add{end+1}=@(theta) theta(j+1)/tau;
  %    models(i).labels=[models(i).labels 6+j];
    end
  end
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
