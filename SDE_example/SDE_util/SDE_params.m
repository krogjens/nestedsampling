function [params] = SDE_params(theta,M) 
n = 2; % Number of paramaters that are always active;
params=zeros(length(M)+n,1);
params(1:n)=theta(1:n);
for j=1:length(M)
  if M(j)==1; %Only transform variable parameters
    n=n+1;
    params(j+2)=theta(n);
  end
end

