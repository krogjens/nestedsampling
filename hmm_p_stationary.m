function pstate_ini = hmm_p_stationary(W,eps,maxi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Calculates the stationary distribution for the
% transition matrix W.
% eps controls the allowed error
% maxi denotes the maximum number of iterations
% pstate_ini is the output column vector
%
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pstate_ini=W(:,1);
A=W*W;
i=0;
while sum(abs(A(:,1)-pstate_ini))>eps && i<maxi
  pstate_ini=A(:,1);
  A=A*A;
  A=A./sum(A,1);
  i=i+1;
end
pstate_ini=A(:,1);
