function [walker_new,delta_a,delta_r,dim_mod]=ns_join_walker_evolve(obs,model,logLstar,walker,delta_u,u_unfix,model1,walker_step)

lenfix=length(u_unfix)-sum(u_unfix);
if sum(u_unfix)>0
  lenobs=(length(delta_u)-length(u_unfix))/sum(u_unfix)+1;
else
  lenobs=1;
end

delta_a=zeros(size(walker.u));;
delta_r=zeros(size(walker.u));;
dim_mod=zeros(size(walker.u));;

walker_new=walker;
walker_try.u=walker.u;

for n=1:lenfix
  walker_try.u(n) = mod(walker_try.u(n) + delta_u(n),1);
end

walker_try.theta=model.invprior(walker_try.u);
walker_try.logl=model.logl(obs,walker_try.theta);

if(walker_try.logl > logLstar)  % Updates if likelihood increased
  delta_a(1:lenfix)=abs(delta_u(1:lenfix));
  walker_new=walker_try;        % Updates walker_new
else
  delta_r(1:lenfix)=abs(delta_u(1:lenfix));
end
dim_mod(1:lenfix)=lenfix;

if lenobs>1 || sum(u_unfix)>0
  walker1.u=NaN(size(u_unfix));
  delta_u1=NaN(size(u_unfix));
  walker1.u(1:lenfix)=walker_new.u(1:lenfix);
  delta_u1(1:lenfix)=zeros(size(1:lenfix));
  index_unfix=(lenfix+1):(lenfix+sum(u_unfix));
  for n=1:lenobs
    index_cur=(lenfix+(n-1)*sum(u_unfix)+1):(lenfix+n*sum(u_unfix));
    walker1.u(index_unfix) = walker_new.u(index_cur);
    delta_u1(index_unfix)=delta_u(index_cur);
    walker1.theta=model1.invprior(walker1.u);
    walker1.logl=model1.logl(obs{n},walker1.theta);
    logLstar1=logLstar-walker_new.logl+walker1.logl;
    [walker1_new,delta_a1,delta_r1,dim_mod1]=walker_step(obs{n},model1,logLstar1,walker1,delta_u1);
    walker_new.u(index_cur)=walker1_new.u(index_unfix);
    delta_a(index_cur)=delta_a1(index_unfix);
    delta_r(index_cur)=delta_r1(index_unfix);
    dim_mod(index_cur)=sum(u_unfix);
    walker_new.logl=walker_new.logl-walker1.logl+walker1_new.logl;
  end
  walker_new.theta=model.invprior(walker_new.u);
  walker_new.logl=model.logl(obs,walker_new.theta);
end
