function SIMPLE_test()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This skeleton routine carries out Bayesian inference with a test of
% convergence for the MCMC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

% There are one model:
% Model 1 is Brownian motion with a drift

%Specify a Jeffreys prior on the diffusion constant
Dmin=10^(-1);
Dmax=10^2;   
%Additionally specify a uniform prior on the drift
mumin=-10000;
mumax=10000;
model.invprior=@(u) [Dmin*exp(u(1)*log(Dmax/Dmin)); (mumax-mumin)*u(2)+mumin];

%Specify options
options.nwalkers=200;   % Number of walkers to be generated
options.stoprat=10^(-1);% Ratio to stop sampling
options.nsteps=30;      % Attempted number of steps in parameter space
%options.ntest=500;
model.options=options;

%Specify the u-generator
model.genu=@() rand(1,2);

%Specify the functions that evolves the walker
%model.evolver=@ns_evolve_rectangle;
%model.evolver=@ns_evolve_exp;

%Specify the logl
log_normal=@(obs,mu,var) -log(sqrt(2*pi*var))*(length(obs)-1)-sum((obs(2:end)-obs(1:(end-1))-mu).^2)/(2*var);
model.logl=@(obs,theta) log_normal(obs,theta(2),2*theta(1));

%Generate the data
data=cumsum(sqrt(2*1)*randn(1000,1));

%Specify the test routine
ntesters=5*options.nwalkers;
model.test=@(obs,model,logLstar,walkers,step_mod) ns_test_evolve_neighbor_min(obs,model,logLstar,walkers,step_mod,ntesters);

%Add defaults
[model,misc]=ns_default_settings(data,model,struct);

%Run the nested sampling algorithm
results = ns_algorithm(data,model);

testlist = results.testlist;
%Plot the result of the test
ncurves=length(testlist);
curves=zeros(length(testlist(1).res)+1,ncurves);
for i=1:ncurves
  curves(2:end,i)=testlist(i).res;
end
meanend=mean(curves(ceil(end/3):end,:));
for i=1:ncurves
  curves(:,i)=2*curves(:,i)/meanend(i)+i-1;
end
figure(1)
plot(transpose(0:(length(curves(:,1))-1)),curves)
hold on
plot(transpose(0:(length(curves(:,1))-1)),ones(length(curves(:,1)),1)*(2:(ncurves+1)),':')
hold off

