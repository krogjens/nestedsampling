function [pvals]=ns_pvalues(obs,rep,scalar)
% This function calculates p-values for the quantitity scalar using the observations obs and replicated trajectories rep. Scalar can actually be a vector or matrix of scalars too.

obs_scalar=scalar(obs,rep(1).theta);
rep_scalar=scalar(rep(1).obs,rep(1).theta);
pvals=(rep_scalar>obs_scalar);
for i=2:length(rep)
  obs_scalar=scalar(obs,rep(i).theta);
  rep_scalar=scalar(rep(i).obs,rep(i).theta);
  pvals=pvals+(rep_scalar>obs_scalar);
end
pvals=pvals/length(rep);


