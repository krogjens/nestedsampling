function logl = fbm_logl(obs,theta,M)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the likelihood of the sample point (theta) for a given model (M) 
% given the data (obs).
%
% The arguments of the function are
% 
% obs - a 2xT matrix of observations
% theta - a 1x5 vector of parameters signifying the initial point of the
%   MCMC process.
% M - a 1x4 vector specifying the variable parameters of the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T=length(obs);
sigma_d=theta(1);                           %Diffusion constant
mu_emit=[theta(2) ; theta(3)]*ones(1,T);    %Mean drift per frame
sigma_n=theta(4);                           %Noise parameter
H=theta(5);                                 %Hurst Exponent

logl=0; 

if M(4)==0; %If the Hurst parameter, H=0.5. Calculate the likelihood normally
	%Preallocate space
    F=zeros(2,T);
	var=zeros(1,T);
    %Calculate initial exponent and variance
	F(:,1)=obs(:,1)-mu_emit(:,1);
	var(1)=2*sigma_n^2 + sigma_d^2;
    %Calculate likelihood of first step
	logl = logl - sum(F(:,1).^2)/(2*var(1)) - log(2*pi*var(1));
	for t=2:T;
        %Update variance and exponent
		var(t)=var(1)-sigma_n^4/var(t-1);
		F(:,t)=obs(:,t)-mu_emit(:,t)+sigma_n^2/var(t-1)*F(:,t-1);
        %Update likelihood to incorporate current step
		logl = logl - sum(F(:,t).^2)/(2*var(t)) - log(2*pi*var(t));
	end

else %If the model has H!=0.5, calculate the likelihood via the 
     % Durbin-Levinson algorithm
	 
     %Preallocate space for arrays
	gamma=zeros(1,T); 
	sigma2=zeros(1,T);
	phi_last=zeros(1,T);
	phi_next=zeros(1,T);
	mu_cond=zeros(2,T);
	
	%Calculate correlations
	gamma(1)=1/2*sigma_d^2*((1+1).^(2*H) -2*1.^(2*H)+(1-1).^(2*H))-sigma_n^2;
	k=2:T;
	gamma(k)=1/2*sigma_d^2*((k+1).^(2*H) -2*k.^(2*H)+(k-1).^(2*H));
	
	%Specify initial variation and similarity
	sigma2(1)=sigma_d^2+2*sigma_n^2;
	phi_last(1)=gamma(1)/sigma2(1);
	mu_cond(:,1)=mu_emit(:,1);
	logl = - sum((obs(:,1)-mu_cond(:,1)).^2)/(2*sigma2(1)) - log(2*pi*sigma2(1));

	%Recursive determination of conditional variance and mean
	for t=2:T;
		%Update variance
		sigma2(t)=sigma2(t-1)*(1-phi_last(t-1)^2);				 
		%Calculate conditional means                                     
		mu_cond(1,t)=mu_emit(1,t)+sum(phi_last(1:t-1).*(obs(1,t-(1:t-1))-mu_emit(1,t-1)));
		mu_cond(2,t)=mu_emit(2,t)+sum(phi_last(1:t-1).*(obs(2,t-(1:t-1))-mu_emit(2,t-1)));
		%Calculate contribution to likelihood	
		logl = logl - sum((obs(:,t)-mu_cond(:,t)).^2)/(2*sigma2(t)) - log(2*pi*sigma2(t));
		%Calculate new similarity
		phi_next(t)=(gamma(t)-sum(phi_last(1:t-1).*gamma(t-(1:t-1))))/sigma2(t); 
	        phi_next(1:t-1)=phi_last(1:t-1)-phi_last(t-(1:t-1))*phi_next(t);	 
		%Save similarity for next iteration
		phi_last=phi_next;		
	end		

end
