function out = hmm_pbm(obs,mu,sigma,out_wish)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function that calculates probability, Rosenblatt transformation
% and its inverse for a single step in each state for a Hidden Markov model
% with Brownian motion in each state.
% obs is a 1xDim vector (Dim=dimensionality of space)
% mu is a MxDim matrix (M=#states)
% sigma is a Mx1 vector
% out_wish is either 'l' (loglikelihood), 'u' (u-values), 'x' (synthetic obs)
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% allow a quick way to have zero mean
if numel(mu)==0
  mu=zeros(size(sigma,1),size(obs,2));
end

if out_wish=='l'
  var=sigma.^2;
  stepdevs=obs-mu;
  out=-log(2*pi*var)*size(stepdevs,2)/2-sum(stepdevs.^2,2)./(2*var);
elseif out_wish=='u'
  stepdevs=obs-mu;
  out=1-erfc(stepdevs./(sqrt(2)*sigma))/2;
elseif out_wish=='x'
  stepdevs = sqrt(2)*sigma.*erfcinv(2-2*obs);
  out = stepdevs+mu;
end
end
