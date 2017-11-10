function [walker_new,delta_a,delta_r]=ns_walker_step(obs,model,logLstar,walker,delta_u)

for n=1:length(walker.u)
  walker_new.u(n) = mod(walker.u(n) + delta_u(n),1);
end
walker_new.theta=model.invprior(walker_new.u);
walker_new.logl=model.logl(obs,walker_new.theta); % Calculates new likelihood

if(walker_new.logl > logLstar)              % Updates if likelihood increased
  delta_a = abs(delta_u);
  delta_r = zeros(size(delta_u));        
else
  walker_new=walker;
  delta_a = zeros(size(delta_u));        
  delta_r = abs(delta_u);
end

