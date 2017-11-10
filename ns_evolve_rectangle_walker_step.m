function [walker_new,step_mod]=ns_evolve_rectangle_walker_step(obs,model,logLstar,walker,step_mod,walker_step)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uses the MCMC technique to find a new sample uniformly distributed inside
% the region of parameter space with higher likelihood than the minimum 
% requirement (logLstar). 
%
% Some of the arguments of the function are
% 
% obs - a 2xT matrix of observations
% walker - the walker that constitutes the starting point of the
%   MCMC process.
% logLstar - the minimum requirement for the likelihood of the new sample
% step_mod - a scalar value that regulates the average step length of the
%   subsequent call of the function. When the remaning parameter space
%   becomes small, the MCMC steps are adjusted in length to ensure a
%   success rate of the MCMC steps of about 50%.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize step_mode if run for the first time
if step_mod==0
  step_mod=ones(size(walker.u));
end

%Specify likelihood function
logl = model.logl;

% Attempt to generate new walker through independent guess
walker_new.u = model.genu();
walker_new.theta = model.invprior(walker_new.u);
walker_new.logl=logl(obs,walker_new.theta); % Calculates new likelihood 

if(walker_new.logl <= logLstar)	% Do MCMC if likelihood requirement failed
   % Find a new walker via a MCMC process.
   % Return input values and logLstar if no steps succeed
   walker_new = walker;
   walker_new.logl=logLstar;

   i = 0;
   reject = zeros(size(walker.u));
   accept = zeros(size(walker.u));
   for i=1:model.options.nsteps
      delta_u=(rand(size(walker.u)) - 0.5) .* step_mod; % Propose step for parameters
      [walker_new,delta_a,delta_r]=walker_step(obs,model,logLstar,walker_new,delta_u);
      accept = accept + delta_a;
      reject = reject + delta_r;
   end
   frac=0.5;
   for n=1:length(step_mod)
     step_mod(n)= 0.5*(1+((1/frac*accept(n)+step_mod(n)/4)/(accept(n)+reject(n)+step_mod(n)))^(1/length(step_mod)))*step_mod(n);
     step_mod(n) = min(step_mod(n),1);	% Update step modifier
   end
end
