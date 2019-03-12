function [percentiles,param_mean,param_stddev,maxLpar] = ns_analyze(samples,model,misc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates percentiles, maximum likelihood parameters,
% mean and standard deviation for the parameters
% in samples. If model.add contains a function of the parameters, then
% percentiles etc. is also calculated for them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   n_samp = length(samples);
   if isfield(model,'ntheta') && length(model.ntheta)>0
     n_theta=model.ntheta;
   else
     n_theta = length(samples(1).theta);
   end
   thetas=NaN(n_theta,n_samp);
   for i=1:n_theta
     for j=1:n_samp
       thetas(i,j)=samples(j).theta(i);
     end
   end
   if isfield(model,'add')
     for j=1:length(model.add)
       add=NaN(1,n_samp);
       for k=1:n_samp
         add(k)=model.add{j}(samples(k).theta);
       end
       thetas=vertcat(thetas,add);
     end
   end
   n_theta=length(thetas(:,1));
   posterior=[samples(1:n_samp).post];
   p_mean = zeros(1,n_theta);
   p_var = zeros(1,n_theta);
   for j=1:n_theta;
     p_mean(j)=sum(thetas(j,:).*posterior);
     p_var(j)=sum((thetas(j,:)-p_mean(j)).^2.*posterior);
   end
   param_mean = p_mean(:);
   param_stddev = sqrt(p_var(:));
   percentiles=ns_percentiles(misc.percentiles_at,thetas,posterior);
   maxLpar=thetas(:,end);
