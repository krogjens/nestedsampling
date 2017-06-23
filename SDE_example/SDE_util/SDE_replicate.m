function obs = SDE_maketrack(obs,theta,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We assume that
% dx/dt = f*mob + kTd/dx*mob + sqrt(2kT*mob)\xi(t)
% with
% mob = D*abs(x)^a
% Time is rescaled; dt -> n*dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = obs(1);
T = length(obs);
D = theta(1);  
a = theta(2); 
f = theta(3);    
var_n = 0; %Localization error (^2) (loc. error)

obs(1) = x0;                                       %Starting point
mob = D * abs(obs(1))^a;                           %Initial mobility
mean_x = f * mob + sign(obs(1)) * a * D * abs(obs(1))^(a-1);      %Initial mean increment
mean_x = mean_x * n;                               %Scale mean increment by scale factor
var_d = 2 * mob;                                   %Initial variance
var_tot = 2 * var_n + var_d * n;                   %Include loc. error
obs(2) = obs(1) + mean_x + sqrt(var_tot) * randn;  %Generate next point
F = obs(2) - obs(1) - mean_x;                      %Calculate deviation
logl = -log(sqrt(2*pi*var_tot)) - F^2/(2 * var_tot);

for t=3:T
   mob = D * abs(obs(t-1))^a;                      %New mobility
   var_d = 2 * mob;                                %New variance
   mean_x = f * mob + sign(obs(t-1)) * a * D * abs(obs(t-1))^(a-1); %New mean
   mean_x = mean_x * n;                            %Scale mean increment by scale factor
   mean_x = mean_x - var_n/var_tot * F;            %Include loc. error through former deviation
   var_tot = var_d * n + 2*var_n - var_n^2/var_tot;    %New variance
   obs(t) = obs(t-1) + mean_x + sqrt(var_tot) * randn;   %Include loc. error
   F = obs(t) - obs(t-1) - mean_x;                 %Calculate deviation
   logl = logl - log(sqrt(2*pi*var_tot)) - F^2/(2 * var_tot);
end
clean = obs;
var_n = theta(4);
obs = obs + randn(1,length(obs)) * sqrt(var_n);
