function [theta] = HMM_invprior(u,ranges,nstates,ndrifts)
%we use a uniform prior for diffusion and drift
%uni=@(u,theta_min,theta_max) (theta_max - theta_min)*u +theta_min;
%and and exponential prior for with mean 1/lambda for transition matrix
%indices
%exp=@(u,lambda) (-log(1-A)*lambda)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The ordering of the parameters in the list will be 
% 1. n diffusion parameters
% 2. j drift parameters
% 3. state transition parameters, which are indices of a transition matrix 
% BUt only the off diagonals are active variables
% example of 3-state, 2-drift hmm
%thetas = [sigma1,sigma2,sigma3,mu1x,mu1y,mu2x,mu2y,T_11,T_21,T_31,T_12,T_22,T_32,T_13,T_23,T_33, P_1,P_2,P_3]
%By passing along the number of states, and drift containing states, we can
%manage all the parameters as a list.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First nstates variables are Diffusion constants (or B) 
%uniform prior on each constrained as B_1 < B_2 < ... B_J for all J states
%otherwise there will be no ordering of states, and multiple identical
%models will count toward the evidence. 
%We might need to forgo that ordering
theta=zeros(1,length(u));
%The first nstates indices are diffusion constants
%with a uniform prior
theta(1:nstates)=(ranges(1,2) - ranges(1,1))*u(1:nstates) +ranges(1,1);

%next 2*ndrift parameters are drift parameters
%also a uniform prior
if ndrifts>0
    theta(nstates+1:nstates+2*ndrifts)=(ranges(2,2) - ranges(2,1))*u(nstates+1:nstates+2*ndrifts) +ranges(2,1);
end
% The final variables will be the transition rates 
%Here we use an exponential prior distribution with mean mu_exp because we are only
%interested in longer lived states.

%first we must reconstruct the transition matrix A by reshaping the u
%values
if nstates >1
    Au=reshape(u(nstates+1+2*ndrifts:end-nstates),nstates,nstates);

    mu_exp=1; %the mean value of the exponential prior, can be adjusted as needed

    %We convert the generated u-values to parameters space via an inverse exponential
    %distribution
    AParamSpace=-log(1-Au)*mu_exp;
    %We adjust the diagonals of the transition matrix 
    AParamSpace=AParamSpace-repmat(sum(AParamSpace,2),1,nstates).*eye(nstates);
    %we then respace for instertion in theta again
    theta(nstates+1+2*ndrifts:end-nstates)=reshape(AParamSpace,1,nstates*nstates);
    
    %finaly the starting probabilities, also uniform prior 
    % to obtain a uniform distribution on an n-dimensional simplex we need
    % a dirichlet distribution, which reduces to an exponential in the
    % uniform case
    
    p = -log(u(end-nstates+1:end));
    S = sum(p);
    p = p/S;
    %I dont think above can be performed by only selecting nstates-1 random
    %values and computing the last one
    
    theta(end-nstates+1:end)=p;
end
