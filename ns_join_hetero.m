function modelsjoin = ns_join_hetero(model,u_unfix,theta_unfix,add_unfix,lenobs)
% NB: ns_join expects this to run if lenobs=1 and all unfix are zeros even though there are multiple trajectories.

modelsjoin = model;

lengthu=length(u_unfix)+(lenobs-1)*sum(u_unfix);
modelsjoin.genu = @() util_generate_u(lengthu);

modelsjoin.logl = @(obs,theta) sum(arrayfun(@(i) model.logl(obs{i},ns_join_reduce(theta,theta_unfix,i)),1:length(obs)));

modelsjoin.invprior =@(u) ns_join_invprior(model.invprior,u,u_unfix,theta_unfix);

for i=1:lenobs
  modelsjoin.labels(ns_join_index(theta_unfix,i)) = model.labels(1:length(theta_unfix));
end

if isfield(model,'scaling')
  modelsjoin.scaling = @(obs,n) arrayfun(@(i) model.scaling(obs{i},n),1:length(obs),'UniformOutput',false);
end
if isfield(model,'replicate')
  modelsjoin.replicate = @(obs,theta,n) arrayfun(@(i) model.replicate(obs{i},ns_join_reduce(theta,theta_unfix,i),n),1:length(obs),'UniformOutput',false);
end
if isfield(model,'logl_n')
  modelsjoin.logl_n = @(obs,theta,n) sum(arrayfun(@(i) model.logl_n(obs{i},ns_join_reduce(theta,theta_unfix,i),n),1:length(obs)));
end
if isfield(model,'add')
  modelsjoin.add={};
  for i=1:length(add_unfix)
    if add_unfix(i)==1
      for j=1:lenobs
        modelsjoin.add{end+1} = @(theta) model.add{i}(ns_join_reduce(theta,theta_unfix,j));
        modelsjoin.labels=[modelsjoin.labels model.labels(length(theta_unfix)+i)];
      end
    else
      modelsjoin.add{end+1} = @(theta) model.add{i}(ns_join_reduce(theta,theta_unfix,1));
      modelsjoin.labels=[modelsjoin.labels model.labels(length(theta_unfix)+i)];
    end
  end 
end
if isfield(model,'checks')
%  modelsjoin.checks = model.checks;
  for i=1:length(model.checks)
    if isfield(model.checks(i).misc,'join')
      modelsjoin.checks(i).scalar =  @(obs,theta) model.checks(i).misc.join(arrayfun(@(j) model.checks(i).scalar(obs{j},ns_join_reduce(theta,theta_unfix,j)),1:length(obs),'UniformOutput',false));
      clear model.checks(i).misc.join;
    else
      if isfield(model.checks(i).misc,'rows') || lenobs==1
        modelsjoin.checks(i).scalar =  @(obs,theta) ns_join_sum(arrayfun(@(j) model.checks(i).scalar(obs{j},ns_join_reduce(theta,theta_unfix,j)),1:length(obs),'UniformOutput',false));
      else
        if isfield(model.checks(i).misc,'columns')
          modelsjoin.checks(i).misc.rows=1:lenobs;
          modelsjoin.checks(i).scalar = @(obs,theta) cell2mat(transpose(arrayfun(@(j) model.checks(i).scalar(obs{j},ns_join_reduce(theta,theta_unfix,j)),1:length(obs),'UniformOutput',false)));
        else
          modelsjoin.checks(i).misc.columns=1:lenobs;
          modelsjoin.checks(i).scalar = @(obs,theta) cell2mat(arrayfun(@(j) model.checks(i).scalar(obs{j},ns_join_reduce(theta,theta_unfix,j)),1:length(obs),'UniformOutput',false));
        end
      end
    end
  end
end

