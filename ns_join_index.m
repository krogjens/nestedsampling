function i_list = ns_join_index(unfix,i)

lenfix=length(unfix)-sum(unfix);
i_list=NaN(size(unfix));
i_list(~logical(unfix)) = 1:lenfix;
i_list(logical(unfix)) = (lenfix+(i-1)*sum(unfix)+1):(lenfix+i*sum(unfix));

