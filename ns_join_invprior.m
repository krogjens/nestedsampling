function theta = ns_join_invprior(invprior,u,u_unfix,theta_unfix)

if sum(u_unfix)>0
  lenobs=(length(u)-length(u_unfix))/sum(u_unfix)+1;
else
  lenobs=1;
end
length_theta = length(theta_unfix)+(lenobs-1)*sum(theta_unfix);
theta=NaN(length_theta,1);

for i=1:lenobs
  theta(ns_join_index(theta_unfix,i)) = invprior(ns_join_reduce(u,u_unfix,i));
end

