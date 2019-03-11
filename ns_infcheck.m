function prob = ns_infcheck(obs,rep,model)

options = model.options;
n_list = options.nlist; 
logl = model.logl;
if isfield(model,'scaling')
   scaling = model.scaling;
   logl_n = model.logl_n;
else
   scaling = @(x,n) x;
end

for n=1:length(n_list)
    n2=n_list(n);
    if n2~= 1 && ~isfield(model,'scaling')
        fprintf('Error: No data rescaling function specified.\n');
    end
    scaled_obs = scaling(obs,n2); % "Scale" original data
    for i=1:length(rep)
        theta = rep(i).theta;
        new_obs = scaling(rep(i).obs,n2);         
        if n2 == 1
            logl_rep = logl(new_obs,theta);
            logl_ref = logl(scaled_obs,theta);
        else
            if ~isfield(model,'logl_n')
               fprintf('Error: No likelihood function specified for rescaled data.\n');
            end
            logl_rep = logl_n(new_obs,theta,n2); 
            logl_ref = logl_n(scaled_obs,theta,n2);
        end 
        H_star(i,n) = logl_rep; 
        H_ref(i,n) = logl_ref;
    end
    prob(n) = mean(H_star(:,n) > H_ref(:,n));
end
