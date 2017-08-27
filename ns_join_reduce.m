function uu = ns_join_reduce(u,unfix,i)

lenfix=length(unfix)-sum(unfix);
unfix_i=find(unfix);
fix_i=find(ones(size(unfix))-unfix);
uu=NaN(size(unfix));

uu(fix_i)=u(1:lenfix);
i_list = (lenfix+(i-1)*sum(unfix)+1):(lenfix+i*sum(unfix));
uu(unfix_i)=u(i_list);

