function modelsjoin = ns_join_simple(model,lenobs)
% If lenobs=1 then it is assumed that all parameters should be common for all trajectories
% If lenobs>1 then no parameters are common for the trajectories

if lenobs==1
  u_unfix=zeros(size(model.genu()));
  theta_unfix=zeros(size(model.invprior(u_unfix+0.5)));
  if isfield(model,'add')
    add_unfix=zeros(1,length(model.add));
  else
    add_unfix=[];
  end
else
  u_unfix=ones(size(model.genu()));
  theta_unfix=ones(size(model.invprior(u_unfix*0.5)));
  if isfield(model,'add')
    add_unfix=ones(1,length(model.add));
  else
    add_unfix=[];
  end
end

modelsjoin = ns_join_hetero(model,u_unfix,theta_unfix,add_unfix,lenobs);

