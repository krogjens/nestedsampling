function logl = SDE_logl_m(obs,theta,n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We assume that
% dx/dt = f*mob + kTd/dx*mob + sqrt(2kT*mob)\xi(t)
% with
% mob = D*abs(x)^a
% Time is rescaled; dt -> n*dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = length(obs);  
D = theta(1);    
a = theta(2);    
f = theta(3);    
var_n = theta(4); %Localization error (^2)

logl = 0; 
mob = D * abs(obs(1))^a;
mean_x = f * mob + sign(obs(1)) * a * D * abs(obs(1))^(a-1);
mean_x = mean_x * n;
F = obs(2) - obs(1) - mean_x;
var_d = 2 * mob;
var_tot = 2 * var_n + var_d * n;
logl = logl - log(sqrt(2*pi*var_tot)) - F^2/(2 * var_tot);

for t=3:T
   mob = D * abs(obs(t-1))^a;
   var_d = 2 * mob;
   mean_x = f * mob + sign(obs(t-1)) * a * D * abs(obs(t-1))^(a-1);
   mean_x = mean_x * n;
   mean_x = mean_x - var_n/var_tot*F;
   F = obs(t) - obs(t-1) - mean_x;
   var_tot = var_d * n + 2*var_n - var_n^2/var_tot;
   logl = logl - log(sqrt(2*pi*var_tot)) - F^2 /(2 * var_tot);
end

