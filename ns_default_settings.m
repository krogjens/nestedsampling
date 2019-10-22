function [models,misc] = ns_default_settings(obs,models,misc)

if isfield(misc,'nssummary')
  if ~isfield(misc,'percentiles_at')
    misc.percentiles_at=[0.02 0.16 0.5 0.84 0.98];
  end
  if ~isfield(misc,'data_id')
    misc.data_id='';
  end
  if ~isfield(models(1),'labels')
    for i=1:length(models)
      models(i).labels=1:length(models(i).genu())
    end
  end
end

for i=1:length(models)
  if ~isfield(models(i),'options')
    models(i).options={};
  end
%  nparams=length(models(i).genu());
  if ~isfield(models(i).options,'nwalkers')
    models(i).options.nwalkers=200; % Number of walkers to be generated
  end
  if ~isfield(models(i).options,'stoprat')
    models(i).options.stoprat=10^(-3); % Z-ratio below which to stop nested sampling
  end
  if ~isfield(models(i).options,'nsteps')
    models(i).options.nsteps=30; % Number of steps per MCMC update
  end
  if ~isfield(models(i).options,'ntest')
    models(i).options.ntest=500;
  end
  if ~isfield(models(i),'evolver') || length(models(i).evolver)==0
    models(i).evolver=@(obs,model,logLstar,walker,step_mod)ns_evolve_separately(obs,model,logLstar,walker,step_mod);
%    models(i).evolver=@(obs,model,logLstar,walker,step_mod)ns_evolve_rectangle(obs,model,logLstar,walker,step_mod);
  end
end

