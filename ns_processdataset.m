function [results]=ns_processdataset(obs,models,misc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs the nested sampling algorithm on the observed data
% in 'obs' on the 'models' with some options stored in 'misc'
% and returning 'results'. More specifically we have:
%
% models - a list of model-structs each with fields
%   logl - the log-likelihood function logl(obs,theta)
%     where theta can be assigned as theta=invprior(walker.u).
%   invprior - the inverse of the cumulative prior;
%     it takes a 1 x options.lengthu vector and returns another vector.
%   options - a struct with at least the fields:
%     lengthu - number of independent parameters in the model
%     nwalkers - number of walkers
%     stoprat - the ratio of evidence in the remaining walkers
%       versus previously summed evidence at which nested sampling stops
%     nsteps - number of attempted steps in a run of ns_evolve
%   labels - a list with the names of the parameters
%   add (optional) - a cell array with functions of theta
% 
% misc - a struct with fields
%   data_id - the first part of the filenames for data and output
%   percentiles_at - a list of values to find percentiles at
%   labels - a list with names for the parameters (theta)
%   titles - further titles/names for ns_print
%
% results - a list of structs with fields
%   logZ - the log of the evidence
%   H - the information
%   Z_norm - the posterior probability of the model
%   logZ_error - estimated error on logZ
%   samples - a list of structs with fields
%     theta - invprior(walker.u) for some walker
%     post - the posterior probability of the sample
%     logl - the log-likelihood models(i).logl(obs,sample.theta)
%   param_mean - estimated averages for parameters (theta)
%   param_stddev - ditto deviations
%   percentiles - the percentiles of theta matching percentiles_at
%   maxLpar - maximum likelihood parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic

% Run nested sampling algorithm for each model
parfor i=1:length(models);
    [results(i).logZ,results(i).H,results(i).samples]...
        =ns_algorithm(obs,models(i));
end

% Calculate total evidence for the models
logZ_tot = log(0);
for i=1:length(models);
    logZ_tot = ns_logsumexp2(logZ_tot,results(i).logZ(1));
end

% Calculate normalized evidence for the models and more
for i = 1:length(models);
    results(i).Z_norm = exp(results(i).logZ(1) - logZ_tot);
    results(i).logZ_error = sqrt(results(i).H(1)/models(i).options.nwalkers);
    [results(i).percentiles,results(i).param_mean,results(i).param_stddev,results(i).maxLpar]...
        =ns_analyze(results(i).samples,models(i),misc);
end

% Calculate probability of  H(replicated obs) > H(obs)
if isfield(models,'replicate')
for i = 1:length(models)
   results(i).prob = ns_infcheck(obs,models(i),results(i).logZ,results(i).samples);
end
else
    for i = 1:length(models)
        results(i).prob = 0.5;
    end
end

%Print a summary of the results to a text file if wanted
if isfield(misc,'nssummary')
  ns_print(results,models,misc)
end

disp('Data processing complete');
toc

