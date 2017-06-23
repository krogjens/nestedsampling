function [obs] = SDE_scaling(obs,n)
n_steps = floor((length(obs)-1)/n);
obs = obs(1 + (0:n_steps)*n);
