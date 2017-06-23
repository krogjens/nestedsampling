function prob = ns_infcheck(obs,model,logZ,samples)

options = model.options;
track_max = options.trackmax;
n_max = options.nmax; 
replicate = model.replicate;
logl = model.logl;
if isfield(model,'scaling')
   scaling = model.scaling;
   logl_n = model.logl_n;
end

post = [samples.post];
post_cum = cumsum(post);
for n=1:n_max
    scaled_obs = scaling(obs,n); % "Scale" original data
    for i=1:track_max
        % Use rand to draw random sample and generate a track 
        point = rand;
        draw = find(post_cum > point,1);
        theta = samples(draw).theta;
        new_obs = replicate(scaled_obs,theta,n);  
        if n == 1
            logl_rep = logl(new_obs,theta);
            logl_ref = logl(scaled_obs,theta);
        else
            logl_rep = logl_n(new_obs,theta,n); 
            logl_ref = logl_n(scaled_obs,theta,n);
        end 
        H_star(i,n) = logl_rep - logZ(n); 
        H_ref(i,n) = logl_ref - logZ(n);
    end
    prob(n) = mean(H_star(:,n) > H_ref(:,n));
end
