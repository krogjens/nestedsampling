%hmm(paramstruct,x,y,T,J,l2p,thetas)
function logl = HMM_logl(obs,theta, nstates,ndrifts,T,dt)

% Prepare arguments for C function call, which involves converting
% transition rate matrix into probabilities as A_prop= expm(dt*A_rate)) and
% also taking log(A_prop) so it it calculated outside of a loop, for
% effeciency
if nstates > 1
    p=theta(end-nstates+1:end);
    logA=log(expm(dt*reshape(theta(nstates+1+2*ndrifts:end-nstates),nstates,nstates)));
else 
    p=1;
    logA=1;
end
B=theta(1:nstates);
u=zeros(nstates*2);
if ndrifts > 0 
    u(1:ndrifts*2)=theta(nstates+1:nstates+2*ndrifts);
end
logl=ForwardAlgorithm(obs,p,logA(:),B,u,nstates,ndrifts,T);
% 
% %hmm_logl(p,logA,B,u,D,x,y, T,J, l2p,thetas)
%     log2pi=log(2*pi);% this can in principle be calculated further out of the loop, but you have to then carry it around
%     sigma=theta(1:nstates);
%     
%     dx=obs(:,1);
%     dy=obs(:,2);
% 
%     if ndrifts > 0
%         mu_x = zeros(nstates); 
%         mu_y = zeros(nstates);
%         mu_x(1:ndrifts)=theta(nstates+1:nstates+ndrifts);  %and then it should be 0 for the remaining states 
%         mu_y(1:ndrifts)=theta(nstates+1+ndrifts:nstates+2*ndrifts);
%     else
%         mu_x = zeros(nstates); 
%         mu_y = zeros(nstates);
%     end
%     
%     if nstates==1 %if only 1 state
%         %The emission model is 2d gaussian with or without drift in
%         %log-space
%         logl=sum(VecTwoDGausslogl(dx,dy,sigma,mu_x,mu_y,log2pi));
%     else
%         p=theta(end-nstates+1:end);
%         %forward algorithm 
%         % if multiple states, use forward algorithm to calculate logl
%         %first time step
%         %logA array with likelihood for each state 
%         logA=log(p)+VecTwoDGausslogl(dx,dy,sigma,mu_x,mu_y,log2pi); %correct here 
%  
%         for t=2:T  %for rest of timesteps
%            B_ik=VecTwoDGausslogl(dx,dy,sigma,mu_x,mu_y, log2pi); %precalculate for each state 
%             for j=1:nstates %over all states for 1 timestep
%                 gamma=LOGPLUS(logA(1)+logA(j,1),logA(2)+logA(j,2));
%                 for i=3:nstates %calculating log sum
%                     gamma=LOGPLUS(gamma,logA(i)+logA(j,i));
%                 end
%                 logA(j)=gamma+B_ik(j);
%             end
%         end
%         %now do logsum
%          logl=LOGPLUS(logA(1),logA(2));
%         for j=3:nstates
%             logl=LOGPLUS(logl,logA(j));
%         end
%     end

end