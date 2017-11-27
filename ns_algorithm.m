function [logZ,H,samples,testlist]=ns_algorithm(obs,model)%
%,logl,logl_n,invprior,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determines the evidence for a model with likelihood function logl
% via the Nested Sampling algorithm. Some variables used below are
%
% walkers - a list of walker-structs each with fields
%   walker.u - a u-value
%   walker.theta - equals invprior(walker.u)
%   walker.logl - equals logl(obs,walker.theta)
% step_mod - a variable that regulates the average step lengths of the
%   MCMC walk. When the remaining parameter space
%   becomes small, the MCMC steps are adjusted in length to ensure a
%   decent success rate of the MCMC steps.
%   If step_mod = 0, then the ns_evolve routine should initiate the variable by itself.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(model.options,'maxsamples')
  maxsampm1=model.options.maxsamples-1;
else
  maxsampm1=Inf;
end

testlist={};

options = model.options;
logl = model.logl;

if isfield(model,'logl_n')
   logl_n = model.logl_n;
end
invprior = model.invprior;

if isfield(options,'nlist')
   nlist = options.nlist;
else
   nlist = 1;
end

logZ=log(0.0)*ones(1,length(nlist)); % Initial evidence
H = zeros(1,length(nlist));          % Initial information

samples = [];
conv_res =[];

%Generate the initial set of walkers with finite likelihood
tries = 0;
logsumwalkers=-Inf;
for i=1:options.nwalkers
  walkers(i).logl = -Inf;
  while walkers(i).logl == -Inf
     tries = tries + 1; % Count the number of total tries to find finite logl samples
     walkers(i).u=model.genu();
     walkers(i).theta=invprior(walkers(i).u);
     walkers(i).logl=logl(obs,invprior(walkers(i).u));
  end
  logsumwalkers=ns_logsumexp2(logsumwalkers,walkers(i).logl);
end

%Outermost interval of prior mass, after one shrink, adjusted for samples with zero likelihood
logwidth=-log(options.nwalkers)+log(options.nwalkers / tries);

%Current ratio of "slab" to total integral and value for stopping
Zrat=Inf; 

step_mod = 0; 		%Tell the ns_evolve routine to initialize step_mod 

i = 1;

while (Zrat>options.stoprat) 	%Stops when the increments of the integral are small

	logwidth=logwidth-log(1.0+1.0/options.nwalkers);   %Shrink interval

	%Identify worst likelihood
	[worst_L,worst]=min([walkers(:).logl]); 

	%Calculate weight of worst walker
	logWt=logwidth+worst_L;

	%Store worst walker in samples
    sample.theta=walkers(worst).theta;
    sample.logl=walkers(worst).logl;
    sample.post=logWt;
    if length(samples)<maxsampm1
      samples = [samples sample];
    else
      k=randi(maxsampm1);
      sumpost=ns_logsumexp2(sample.post,samples(k).post);
      if log(rand)<sample.post-sumpost
        samples(k)=sample;
      end
      samples(k).post=sumpost;
    end

	%Update evidence and check for stopping criteria
	logZnew=ns_logsumexp2(logZ(1),logWt); 	% Updates evidence
    if i == 1
        H(1) = exp(logWt - logZnew) * worst_L - logZnew;
    else
        H(1) = exp(logWt - logZnew) * worst_L + exp(logZ(1) - logZnew) * (H(1) + logZ(1)) - logZnew;
    end
    logZ(1) = logZnew;
    if isfield(model,'scaling')
        for n = 2:length(nlist);
            n2 = nlist(n);
            sc_obs = model.scaling(obs,n2);
            worst_Ln = logl_n(sc_obs,invprior(walkers(worst).u),n2);
            logWt = logwidth + worst_Ln;
            logZ(n) = ns_logsumexp2(logZ(n),logWt);
        end
    end  
	
	%Find random walker to initiate generation of new walker
	copy = ceil(options.nwalkers*rand);  	%Choose random number 1<copy<n_prior
	while(copy==worst && options.nwalkers>1) 
		copy = ceil(options.nwalkers*rand);
	end
	logLstar=walkers(worst).logl;           %New likelihood constraint

	%Evolve copied walker within constraint
        [walker_new,step_mod]=model.evolver(obs,model,logLstar,walkers(copy),step_mod);
	walkers(worst)=walker_new;           %Insert new walker
        logsumwalkers=ns_logsumexp2(logsumwalkers,walker_new.logl+log(1-exp(worst_L-walker_new.logl)));
        Zrat=exp(logwidth+logsumwalkers-logZ(1));
        if mod(i,model.options.ntest) == 0
          fprintf('After %i iterations with %i parameter(s), Zrat = %.3g\n',i,length(walker_new.u),Zrat);
          if isfield(model,'test')
            testlist(i/model.options.ntest).res=model.test(obs,model,logLstar,walkers,step_mod);
          end
        end
	i = i + 1;
end

%Add the remaning samples to the evidence estimate and sample output
[~,I]=sort([walkers(:).logl]);
walkers=walkers(I);
for j=1:options.nwalkers
	logWt=logwidth + walkers(j).logl; 
	logZnew=ns_logsumexp2(logZ(1),logWt);
    H(1) = exp(logWt - logZnew) * walkers(j).logl + exp(logZ(1) - logZnew) * (H(1) + logZ(1)) - logZnew;
    logZ(1) = logZnew;
    sample.theta=walkers(j).theta;
    sample.logl=walkers(j).logl;
    sample.post=logWt;
    if length(samples)<maxsampm1 || j==options.nwalkers
      samples = [samples sample];
    else
      k=randi(maxsampm1);
      sumpost=ns_logsumexp2(sample.post,samples(k).post);
      if log(rand)<sample.post-sumpost
        samples(k)=sample;
      end
      samples(k).post=sumpost;
    end
    if isfield(model,'scaling')
        for n = 2:length(nlist);
            n2 = nlist(n);
            sc_obs = model.scaling(obs,n2);
            worst_Ln = logl_n(sc_obs,invprior(walkers(j).u),n2);
            logWt = logwidth + worst_Ln;
            logZ(n) = ns_logsumexp2(logZ(n),logWt);
        end
    end  
end

%Calculate posterior probability of the samples
for j=1:length(samples)
    samples(j).post=exp(samples(j).post-logZ(1));
end

%Calculate evidence for scaled tracks to check vs replications
if isfield(model,'scaling')
    for n = 2:length(nlist)
        n2=nlist(n);
        H(n) = -logZ(n);
        sc_obs = model.scaling(obs,n2);
        for j = 1:length(samples)
            tht = samples(j).theta;
            H(n) = H(n) + logl_n(sc_obs,tht,n2)*samples(j).post;
        end
    end
end
