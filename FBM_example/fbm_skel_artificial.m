function fbm_skel_artificial()

misc.data_id = 'subdiffusive_data_and_results/subdiffusive_track';

data_path = [misc.data_id,'.txt'];

%loading data 
data = load(data_path);
data = transpose(data);
%Subtract points to obtain steps
obs = data(:,2:end) - data(:,1:end-1);

log10l_true = fbm_logl(obs,[20 0 10 10 0.25],[1 1 1 1])/log(10)

%Tell ns_processdataset to write a summary file
misc.nssummary=['_results.txt'];

%Specify framerate and pixelsize
p_size = 1;  %Pixelsize in data in micrometer/pixel
tau = 0.5;     %Time between data point in seconds

%Specify prior ranges (note: micrometer since it is the length scale used in the data-file)
sigmaHmin=10^(0);  %Minimum value for step deviation in micrometer
sigmaHmax=10^(3);      %Maximum value for step deviation in micrometer

vmax=1000;      %Maximum absolute value of the bias per step

noisemin=0;  %Minimum value for noise deviation parameter in micrometers
noisemax=1000;  %Maximum value for noise deviation parameter in micrometers

Hmin=0;      %Minimum value for Hurst exponent
Hmax=1;      %Maximum value for Hurst exponent

%Convert prior ranges to dimensions of pixels and frames for convenience 
% ranges - a 2x5 array of minimum and maximum values for the 5 parameter:
%   - diffusion coefficient, drift velocity (in two dimentions), noise
%     parameter, and Hurst exponent.
ranges=[sigmaHmin sigmaHmax;...
    -vmax vmax ;-vmax vmax;...
    noisemin noisemax; Hmin Hmax];

%Specify options
options.nlist= [1 2 4 8 16 32];
options.trackmax = 100;

%Specify the models
final=@(x) x(end); % extract last element in a vector
MM = [0 0 0 0; 1 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 0; 1 1 0 1; 0 0 1 1; 1 1 1 1];
i=8;
for k=1:5;
  nparams=sum(MM(i,:))+1;
  if k==1
    models(k).options.nsteps=30; % Number of steps per MCMC update
    models(k).evolver=@(obs,model,logLstar,walker,step_mod)ns_evolve_separately(obs,model,logLstar,walker,step_mod);
    options.Nparfor = 40;
  elseif k==2
    models(k).options.nsteps=30; % Number of steps per MCMC update
    models(k).evolver=@(obs,model,logLstar,walker,step_mod)ns_evolve_separately(obs,model,logLstar,walker,step_mod);
    options.Nparfor = 1;
  elseif k==3
    models(k).options.nsteps=30*nparams; % Number of steps per MCMC update
    models(k).evolver=@(obs,model,logLstar,walker,step_mod)ns_evolve_rectangle(obs,model,logLstar,walker,step_mod);
    options.Nparfor = 40;
  elseif k==4
    models(k).options.nsteps=30*nparams; % Number of steps per MCMC update
    models(k).evolver=@(obs,model,logLstar,walker,step_mod)ns_evolve_rectangle(obs,model,logLstar,walker,step_mod);
    options.Nparfor = 2;
  elseif k==5
    models(k).options.nsteps=30*nparams; % Number of steps per MCMC update
    models(k).evolver=@(obs,model,logLstar,walker,step_mod)ns_evolve_rectangle(obs,model,logLstar,walker,step_mod);
    options.Nparfor = 1;
  end
  models(k).options=options;
  models(k).genu=@() util_generate_u(nparams);
  models(k).logl=@(obs,theta) fbm_logl(obs,fbm_params(theta,MM(i,:)),MM(i,:));
  models(k).invprior=@(u) fbm_invprior(u,ranges,MM(i,:));
  models(k).scaling = @(obs,n) fbm_scaling(obs,n);
  models(k).replicate = @(obs,theta,n) fbm_replicate(obs,fbm_params(theta,MM(i,:)),n);
  models(k).logl_n = @(obs,theta,n) fbm_logl_n(obs,fbm_params(theta,MM(i,:)),MM(i,:),n);
  models(k).labels=[1];
  for j=1:length(MM(i,:));
    if MM(i,j)==1
      models(k).labels=[models(k).labels j+1];
    end
  end
  models(k).add{1}=@(theta) theta(1)^2/(2*tau^(2*final(fbm_params(theta,MM(i,:)))));
  models(k).labels=[models(k).labels 6];
  for j=1:2
    if MM(i,j)==1
      models(k).add{end+1}=@(theta) theta(j+1)/tau;
      models(k).labels=[models(k).labels 6+j];
    end
  end
end

%Labels for the parameters
misc.labels=...
['Step deviation: ';...
 'x-bias/step:    ';...
 'y-bias/step:    ';...
 'Measurem. err.: ';...
 'Hurst exponent: ';...
 'D_H constant:   ';...
 'x-velocity:     ';...
 'y-velocity:     '];

path=[misc.data_id,'_output'];

[results] = ns_processdataset(obs,models,misc);
save(path,'results')

