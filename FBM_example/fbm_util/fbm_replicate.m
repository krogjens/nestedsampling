function [steps] = fbm_replicate(obs,theta,n)
%[steps] = fbm_replicate(theta,n_steps,n)

d=2; % nummber of dimensions
T = length(obs);
 
sigma_d = theta(1);
mu_emit = [theta(2) ; theta(3)]*ones(1,T);
sigma_n = theta(4);
H = theta(5);

gam(1)= 1/2*sigma_d^2 * n^(2 * H) * (2^(2*H)-2*1) - sigma_n^2;
k=2:T;
gam(k)= 1/2*sigma_d^2 * n^(2 * H) * ( (k+1).^(2*H) -2*k.^(2*H) + (k-1).^(2*H) );

var(1)=sigma_d^2 * n^(2 * H) + 2 * sigma_n^2;
phi_last(1)=gam(1)/var(1);
mu_cond(1:d,1) = mu_emit(1:d,1);


steps(:,1)=sqrt(var(1)) * randn(d,1) + mu_cond(:,1);

    for t=2:T;
        var(t) = var(t-1)*(1-phi_last(t-1)^2);
        for i = 1:d
            mu_cond(i,t) = mu_emit(i,t) + sum(phi_last(1:t-1).*(steps(i,t-(1:t-1)) - mu_emit(i,t-1)));
        end
        steps(:,t) = mu_cond(:,t) + sqrt(var(t))*randn(d,1); 
        phi_next(t) = (gam(t) - sum(phi_last(1:t-1).*gam(t-(1:t-1))))/var(t);
        phi_next(1:t-1)=phi_last(1:t-1)-phi_last(t-(1:t-1))*phi_next(t);
        phi_last = phi_next;
    end


