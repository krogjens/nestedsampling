function rep = ns_replicate(obs,model,samples,n_tracks)

rep = []; %Initiate rep structure
replicate = model.replicate; 
post = [samples.post]; 
post_cum = cumsum(post); % Cumsum posterior distribution 
                         % to enable draws
for i =1:n_tracks;
   point = rand;
   draw = find(post_cum > point,1);
   theta = samples(draw).theta;
   new_obs = replicate(obs,theta,1);
   rep(i).theta = theta;
   rep(i).obs = new_obs;
end
