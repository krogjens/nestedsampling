function [vals]=util_stepdist(obs,nsteps,percentiles)
% Calculates estimated percentiles for the cumulative distribution of steps for obs.
% The steps can be of different length as described by nsteps

vals=zeros(length(nsteps),length(percentiles));
len=length(obs);
for i=1:length(nsteps)
  n=nsteps(i);
  steps=obs((n+1):len)-obs(1:(len-n));
  sorted=sort(steps);
  for j=1:length(percentiles)
    vals(i,j)=sorted(ceil(percentiles(j)*(len-n)));
  end
end

