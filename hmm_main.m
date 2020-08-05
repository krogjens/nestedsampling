function out = hmm_main(obs,func,W,pstate_ini,out_wish)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
% function that implements hidden Markov models for calculation of likelihood,
% Rosenblatt transformation and its inverse.
% obs is a NxDim matrix of observations/u-values (N=#times, Dim=dimensionality)
% func(obs(i,:),out_wish) gives a matrix that contains, when out_wish==
%   'l', a column with log-probability of obs(i,:) given the M different states
%   'u', M rows with Dim u-values corresponding to each of the M states
%   'x', M rows with Dim x-values corresponding to each of the M states
% W(i,j) gives probability of next state being i given current state j
% pstate_ini is a Mx1 vector with the initial probabilities (M=#states)
% out_wish is either 'l' (loglikelihood), 'u' (u-values), 'x' (synthetic obs)
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=size(pstate_ini,1);  % M=#states
logWT=log(W');
out=[];
logpstatexprev=log(pstate_ini); % initialize joint probability of state@i and x(1:(i-1),:)
for i=1:size(obs,1)
  obs_cur=obs(i,:);
  if out_wish~='l'
    logpxprev=bmc_logsumexpcol(logpstatexprev); % probability of x(1:(i-1),:)
    pstate_xprev=exp(logpstatexprev-logpxprev); % probability of state@i given x(1:(i-1),:)
    if out_wish=='u'
      out=[out; pstate_xprev'*func(obs_cur,'u')];
    else
      minfunc=@(x) sum((obs_cur-pstate_xprev'*func(x,'u')).^2);
      obs_cur=fminsearch(minfunc,pstate_xprev'*func(obs_cur,'x'));
      out=[out; obs_cur];
    end
  end
  logpstatex=func(obs_cur,'l')+logpstatexprev; % joint probability of state@i and x(1:i,:)
  logpstatexprev=bmc_logsumexpcol(logWT+logpstatex)';  % for next iteration
end
if out_wish=='l'
  out=bmc_logsumexpcol(logpstatex); % probability of x(1:N,:), i.e., the likelihood of the trajectory
end
end
