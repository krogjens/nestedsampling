function [logZ,H,samples]=ns_algorithm(obs,model)%
%,logl,logl_n,invprior,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the evidence for a model with likelihood function logl
% via the Nested Sampling algorithm. Some variables used below are
%
% walkers - a list of walker-structs each with fields
%   walker.u - a u-value, i.e., a 1 x model.options.lengthu vector
%   walker.theta - equals invprior(walker.u)
%   walker.logl - equals logl(obs,walker.theta)
% step_mod - a variable that regulates the average step lengths of the
%   MCMC walk. When the remaining parameter space
%   becomes small, the MCMC steps are adjusted in length to ensure a
%   decent success rate of the MCMC steps.
%   If step_mod = 0, then ns_evolve should initiate the variable by itself.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


options = model.options;
logl = model.logl;
if isfield(model,'logl_n')
   logl_n = model.logl_n;
end
invprior = model.invprior;

if isfield(options,'nmax')
   nmax = options.nmax;
else
   nmax = 1;
end

logZ=log(0.0)*ones(1,nmax); % Initial evidence
H = zeros(1,nmax);          % Initial information

samples = [];

%Generate the initial set of walkers
for i=1:options.nwalkers
  walkers(i).u=rand(1,options.lengthu);
  walkers(i).theta=invprior(walkers(i).u);
end

%Calculate likelihood for the walkers
for i = 1:options.nwalkers;
        walkers(i).logl=logl(obs,invprior(walkers(i).u));
end
%Outermost interval of prior mass
logwidth=-log(options.nwalkers+1);

%Current ratio of "slab" to total integral and value for stopping
Zrat=1;

step_mod = 0; 		%Tell ns_evolve to initialize step_mod 

i = 1;

while (Zrat>options.stoprat); 	%Stops when the increments of the integral are small

	%Identify worst likelihood
	[worst_L,worst]=min([walkers(:).logl]);

	%Calculate weight of worst walker
	logWt=logwidth+worst_L;

	%Store worst walker in samples
        sample.theta=walkers(worst).theta;
        sample.logl=walkers(worst).logl;
        sample.post=logWt;
        samples = [samples sample];

	%Update evidence and check for stopping criteria
	logZnew=ns_logsumexp2(logZ(1),logWt); 	% Updates evidence
        if i == 1
            H(1) = exp(logWt - logZnew) * worst_L - logZnew;
        else
            H(1) = exp(logWt - logZnew) * worst_L + exp(logZ(1) - logZnew) * (H(1) + logZ(1)) - logZnew;
        end
        logZ(1) = logZnew;
	Zrat=exp(log(options.nwalkers)+logWt-logZ(1));  % Measures increment of evidence
        
        if isfield(model,'scaling')
           for n = 2:nmax;
               sc_obs = model.scaling(obs,n);
               worst_L = logl_n(sc_obs,invprior(walkers(worst).u),n);
               logWt = logwidth + worst_L;
               logZnew = ns_logsumexp2(logZ(n),logWt);
               if i == 1
                  H(n) = exp(logWt - logZnew) * worst_L - logZnew;
               else
                  H(n) = exp(logWt - logZnew) * worst_L + exp(logZ(n) - logZnew) * (H(n) + logZ(n)) - logZnew;
               end
               logZ(n) = logZnew;
            end
         end   


	if Zrat< options.stoprat 		% Make sure the maximum weight is not too large
		[best_L,~]=max([walkers(:).logl]);
		Zrat = exp(log(options.nwalkers)+logwidth + best_L - logZ(1));
	end

	%Find random walker to initiate generation of new walker
	copy = ceil(options.nwalkers*rand);  	%Choose random number 1<copy<n_prior
	while(copy==worst && options.nwalkers>1) 
		copy = ceil(options.nwalkers*rand);
	end
	logLstar=walkers(worst).logl;           %New likelihood constraint

	%Evolve copied walker within constraint
	[walker_new,step_mod]=ns_evolve(obs,logl,invprior,logLstar,walkers(copy),step_mod,options);
	walkers(worst)=walker_new;           %Insert new walker
	logwidth=logwidth-log(1.0+1.0/options.nwalkers);   %Shrink interval
    if mod(i,100) == 0
       fprintf('After %i iterations of nested sampling, Zrat =%.4f\n',i,Zrat);
    end

	i = i + 1;
end

%Add the remaning samples to the evidence estimate and sample output
[~,I]=sort([walkers(:).logl]);
walkers=walkers(I);
for j=1:options.nwalkers;
	logWt=logwidth + walkers(j).logl; 
	logZnew=ns_logsumexp2(logZ(1),logWt);
        H(1) = exp(logWt - logZnew) * walkers(j).logl + exp(logZ(1) - logZnew) * (H(1) + logZ(1)) - logZnew;
        logZ(1) = logZnew;
        sample.theta=walkers(j).theta;
        sample.logl=walkers(j).logl;
        sample.post=logWt;
        samples = [samples sample];
        if isfield(model,'scaling')
           for n = 2:nmax
               sc_obs = model.scaling(obs,n);
               worst_L = logl_n(sc_obs,invprior(walkers(j).u),n);
               logWt = logwidth + worst_L;
               logZnew = ns_logsumexp2(logZ(n),logWt);
               H(n) = exp(logWt - logZnew) * worst_L + exp(logZ(n) - logZnew) * (H(n) + logZ(n)) - logZnew;
               logZ(n) = logZnew;
           end
        end
end

%Calculate posterior probability of the samples
for j=1:length(samples)
  samples(j).post=exp(samples(j).post-logZ(1));
end

