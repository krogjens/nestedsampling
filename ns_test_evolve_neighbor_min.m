function [conv_res]=ns_test_evolve_neighbor_min(obs,model,logLstar,walkers,step_mod,ntesters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tests whether the MCMC technique can find a new sample
% uniformly distributed inside
% the region of parameter space with higher likelihood than the minimum 
% requirement (logLstar). 
%
% Some of the arguments of the function are
% 
% obs - the observations
% walkers - the walkers that constitutes the starting point of the
%   MCMC process.
% logLstar - the minimum requirement for the likelihood of the new sample
% step_mod - a variable that regulates the average step length of the
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

%Find the unit along direction m as the average abs-coordinate difference to the nearest neighbor for the walkers
iif = @(c,t,f) c.*(t-f) + f;
n2IJ = @(I,J) iif(I<J,(I-1).*(nwalkers-I/2)+J-I,(J-1).*(nwalkers-J/2)+I-J);
if nwalkers>1
  dists=pdist(walkcoords);
  unit=zeros(1,lenu);
  for j=1:nwalkers
    [not_used,I]=min(dists(n2IJ(j,[1:(j-1) (j+1):nwalkers])));
    I=iif(I<j,I,I+1);
    for m=1:lenu
      unit(m)=unit(m)+abs( mod(walkcoords(I,m)-walkcoords(j,m)+0.5,1)-0.5 );
    end
  end
  unit=unit/nwalkers;
else
  unit=ones(1,lenu);
end

%Find the minimum of the abs-coordinate differences to the nearest neighbor original walker
%as a function of number of steps
conv_res=zeros(1,teststeps);
for n=1:teststeps
  for i=1:ntesters
    [testers(i),not_used]=model.evolver(obs,model,logLstar,testers(i),step_mod);
    [not_used,I] = pdist2(walkcoords,testers(i).u,'euclidean','Smallest',1);
    mindist=Inf;
    for m=1:lenu
      mindist=min(mindist,abs( mod(walkcoords(I,m)+0.5-testers(i).u(m),1)-0.5 )/unit(m));
    end
    conv_res(n)=conv_res(n)+mindist/ntesters;
  end
end
