function [theta] = fbm_invprior(u,ranges,M)
%jeff=@(u,thmin,thmax) thmin*exp(u*log(thmax/thmin));
%uni=@(u,theta_min,theta_max) (theta_max - theta_min)*u +theta_min;
theta=zeros(length(u),1);
theta(1)=util_jeff(u(1),ranges(1,1),ranges(1,2));
n=1;
for j=1:length(M)
  if M(j)==1; %Only transform variable parameters
    n=n+1;
    theta(n)=util_uni(u(n),ranges(j+1,1),ranges(j+1,2));
  end
end

