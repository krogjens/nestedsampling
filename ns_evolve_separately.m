function [walker_new,step_mod]=ns_evolve_separately(obs,model,logLstar,walker,step_mod)

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
   while(i < model.options.nsteps)
     % Propose step for parameters
     % Ensure that new walker.u is between 0 and 1
     for n=1:length(walker.u)
       if isfield(model,'u_evolve')
         delta_u_vec=zeros(size(walker.u));
         delta_u_vec(n)=(rand - 0.5) .* step_mod(n);
         walker.u = model.u_evolve(walker_new.u,delta_u_vec);
       else
         delta_u=(rand - 0.5) .* step_mod(n);
         walker.u(n) = mod(walker_new.u(n) + delta_u,1);
       end
       walker.theta=model.invprior(walker.u);
       walker.logl=logl(obs,walker.theta); % Calculates new likelihood
       if(walker.logl > logLstar)          % Updates if likelihood increased
         walker_new=walker;	  % Updates walker_new
       else
         reject(n) = reject(n) + 1;
         walker=walker_new;	  % Updates walker_new
       end
     end
     i = i + 1;
   end
   for n=1:length(step_mod)
     R = reject(n) / model.options.nsteps;  % Ratio of rejectance
     step_mod(n) = min(step_mod(n) * exp(0.25-R),1);     % Update step modifier
     if R==1
       fprintf('Warning (%i parameters): an iteration occurred with no moves in direction %i\n',length(walker.u),n);
     end
   end
end

