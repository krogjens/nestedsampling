function percentiles = ns_percentiles(p_list,theta_list,posterior)

ntheta=length(theta_list(:,1));
percentiles=NaN(ntheta,length(p_list));
for j=1:ntheta
  [theta_sort, I] = sort(theta_list(j,:));
  post_sort=posterior(I);

  post_cum=[0 cumsum(post_sort) 1];
  theta_sort=[theta_sort(1) theta_sort theta_sort(end)];

  n=2;
  for m=1:length(p_list)
     while p_list(m)>post_cum(n)
       n=n+1;
     end
     percentiles(j,m)= (theta_sort(n)*(p_list(m)-post_cum(n-1))+theta_sort(n-1)*(post_cum(n)-p_list(m)))/(post_cum(n)-post_cum(n-1));
  end
end

