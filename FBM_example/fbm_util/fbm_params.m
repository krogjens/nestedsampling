function [params] = fbm_params(theta,M)
params=zeros(length(M)+1,1);
params(length(M)+1,1)=0.5;
params(1)=theta(1);
n=1;
for j=1:length(M)
  if M(j)==1; %Only transform variable parameters
    n=n+1;
    params(j+1)=theta(n);
  end
end

