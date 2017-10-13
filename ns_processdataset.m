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
%     it takes a u-value vector and returns another vector.
%   genu - function that generates the u-values.
%   options - a struct with at least the fields:
%     nwalkers - number of walkers
%     stoprat - the ratio of evidence in the remaining walkers
%       versus previously summed evidence at which nested sampling stops
%     nsteps - number of attempted steps in a run of ns_evolve
%     nlist - a list specifying the time rescalings to be used for
%       information content checks
%     trackmax - the amount of tracks to replicate for the information content check
%
%   labels - a list with the names of the parameters
%   add (optional) - a cell array with functions of theta
%   replicate - function that generates artificial data by
%     sampling from the posterier probability distribution for the parameters    
%   scaling (optional) - a function that reshapes the data to what would have been observed
%     with lower sampling rates by removing data points.
%   logl_n (optional) - a likelihood function logl_n(obs,theta,n)
%     like "logl" except that takes in a third argument, n, which states 
%     the amount of time step rescaling used for the inserted observations. 
%
%   If the optional functions are omitted, the nlist must consist only of ones.
%   checks (optional) - a list of sctructs with fields
%     scalar - a function that takes a trajectory (e.g. obs) and returns a matrix of real numbers
%     misc.rows (optional) - a list of numbers describing the rows in scalar
%     misc.colums (optional) - a list of number describing the columns in scalar
%     misc.labels a cell array with 1 or 2 elements
%       1st element: a string describing the model check
%       2nd element: (optional) a string describing columns and/or rows
%
% misc - a struct with fields
%   data_id - the first part of the filenames for data and output
%   percentiles_at - a list of values to find percentiles at
%   labels - a list with names for the parameters (theta)
%   titles - further titles/names for ns_print
%   nssummary (optional) - final part of name of a summary file
%   append (optional) - what to write initially if nssummary should be opened in append mode
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
%
% The routine also outputs percentiles for the parameters as well as
% performing an information content model check with the ns_infcheck function.
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

% Replicate data for all models and perform tests
% (Test)  Calculate probability of  H(replicated obs) > H(obs)
% optionally calculates further model checks according to the scalars in the field checks
if isfield(models,'replicate')
   for i = 1:length(models)
      rep = ns_replicate(obs,models(i),results(i).samples,models(i).options.trackmax);
      results(i).prob = ns_infcheck(obs,rep,models(i));
      if isfield(models,'checks')
        for j=1:length(models(i).checks)
          results(i).checks(j).pvals=ns_pvalues(obs,rep,models(i).checks(j).scalar);
        end
      end
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

