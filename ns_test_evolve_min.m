function [conv_res]=ns_test_evolve(obs,model,logLstar,walkers,step_mod,ntesters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests whether the MCMC technique can find a new sample
% uniformly distributed inside
% the region of parameter space with higher likelihood than the minimum 
% requirement (logLstar). 
%
% Some of the arguments of the function are
% 
% obs - a 2xT matrix of observations
% walkers - the walkers that constitutes the starting point of the
%   MCMC process.
% logLstar - the minimum requirement for the likelihood of the new sample
% step_mod - a scalar value that regulates the average step length of the
%   subsequent call of the function. When the remaning parameter space
%   becomes small, the MCMC steps are adjusted in length to ensure a
%   success rate of the MCMC steps of about 50%.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nwalkers=length(walkers);
teststeps=model.options.nsteps;
model.options.nsteps=1;
lenu=length(walkers(1).u);
ind=randperm(nwalkers);
for i=1:ntesters
 testers(i)=walkers(ind(mod(i,nwalkers)+1));
end
walkcoords=zeros(nwalkers,lenu);
for i=1:nwalkers
  walkcoords(i,:)=walkers(i).u;
end 

if nwalkers>1
  unit=zeros(1,lenu);
  for m=1:lenu
    for j=1:nwalkers
      unit(m)=unit(m)+min(abs( mod(walkcoords([1:(j-1) (j+1):nwalkers],m)-walkcoords(j,m)+0.5,1)-0.5 ));
    end
  end
  unit=unit/nwalkers;
else
  unit=ones(1,lenu);
end

conv_res=zeros(1,teststeps);
for n=1:teststeps
  for i=1:ntesters
    [testers(i),not_used]=model.evolver(obs,model,logLstar,testers(i),step_mod);
    mindist=Inf;
    for m=1:lenu
      mindist=min(mindist,min(abs((mod(walkcoords(:,m)+0.5-testers(i).u(m),1)-0.5)))/unit(m));
    end
    conv_res(n)=conv_res(n)+mindist/ntesters;
  end
end
