function [pvals]=ns_pvalues(obs,rep,scalar)
% This function calculates p-values for the quantitity scalar using the observations obs and replicated trajectories rep. Scalar can actually be a vector or matrix of scalars too.

obs_scalar=scalar(obs);
pvals=zeros(size(obs_scalar));;
for i=1:length(rep)
  rep_scalar=scalar(rep(i).obs);
  pvals=pvals+(rep_scalar>obs_scalar);
end
pvals=pvals/length(rep);


