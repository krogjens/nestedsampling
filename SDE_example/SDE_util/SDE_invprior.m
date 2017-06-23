function [theta] = SDE_invprior(u,ranges,M)
%jeff=@(u,thmin,thmax) thmin*exp(u*log(thmax/thmin));
%uni=@(u,theta_min,theta_max) (theta_max - theta_min)*u +theta_min;
theta=zeros(length(u),1);
theta(1)=jeff(u(1),ranges(1,1),ranges(1,2));
theta(2)=uni(u(2),ranges(2,1),ranges(2,2));
n=2; % Number of parameters that are always active
for j=1:length(M)
  if M(j)==1; %Only transform variable parameters
    n = n + 1; %Update number of active parameters
    theta(n)=uni(u(n),ranges(j+2,1),ranges(j+2,2));
  end
end

