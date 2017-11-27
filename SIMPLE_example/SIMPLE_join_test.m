function SIMPLE_join_test()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This skeleton routine carries out Bayesian inference with a test of
% convergence for the MCMC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% There are one model:
% Model 1 and 2 are Brownian motion with a drift

%Specify a Jeffreys prior on the diffusion constant
Dmin=10^(-1);
Dmax=10^2;   
%Additionally specify a uniform prior on the drift
mumin=-100;
mumax=100;
model.invprior=@(u) [Dmin*exp(u(1)*log(Dmax/Dmin)); (mumax-mumin)*u(2)+mumin];

%Specify options
options.nwalkers=200;   % Number of walkers to be generated
options.stoprat=10^(-2);% Ratio to stop sampling
options.nsteps=30;      % Attempted number of steps in parameter space
options.ntest=500;
options.maxsamples=10000; % Restrict the number of samples to use for analysis and saving
model.options=options;

%Specify the u-generator
model.genu=@() rand(1,2);

%Specify the logl
log_normal=@(obs,mu,var) -log(sqrt(2*pi*var))*(length(obs)-1)-sum((obs(2:end)-obs(1:(end-1))-mu).^2)/(2*var);
model.logl=@(obs,theta) log_normal(obs,theta(2),2*theta(1));

%Specify the index for the labels
model.labels=[1 2];

%Percentiles to be calculated
misc.percentiles_at=[0.02 0.16 0.5 0.84 0.98];

%Labels for the parameters and percentile line in the output text-file
misc.labels=...
['Diffusion constant: ';...
 'Drift:              '];
misc.titles=...
{'Percentile at:      '};

%Specify output filename beginnings
misc.data_id = 'simple';

%Generate the data
lenobs=10;
data={};
for n=1:lenobs
  data{n}=cumsum(sqrt(2*1)*randn(200,1));
end

%Tell ns_print to write a summary-file
misc.nssummary=['_results.txt'];

%Specify the test routine
ntesters=5*options.nwalkers;
model.test=@(obs,model,logLstar,walkers,step_mod) ns_test_evolve_min(obs,model,logLstar,walkers,step_mod,ntesters);

u_unfix=[0 1];
theta_unfix=u_unfix;
add_unfix=[];
model_join=ns_join_hetero(model,u_unfix,theta_unfix,add_unfix,lenobs);


%Specify the special functions that evolves the walker for Model 1
walker1_step=@ns_walker_step;
walker_step=@(obs,model_full,logLstar,walker,delta_u) ns_join_walker_evolve(obs,model_full,logLstar,walker,delta_u,u_unfix,model,walker1_step);
model_join.evolver=@(obs,model,logLstar,walker,step_mod) ns_evolve_rectangle_walker_step(obs,model,logLstar,walker,step_mod,walker_step);
models(1)=model_join;

% Switch to standard evolver for Model 2
model_join.evolver=@ns_evolve_rectangle; 
models(2)=model_join;

%Run the nested sampling algorithm for all models and compile results
[results] = ns_processdataset(data,models,misc);

%Save results
path=[misc.data_id,'_output'];
save(path,'results')
save(path,'data','-append')

%Plot the result of the test
figure(1)
ns_test_evolve_plot(results(1).testlist)

%Plot the result of the test
figure(2)
ns_test_evolve_plot(results(2).testlist)

