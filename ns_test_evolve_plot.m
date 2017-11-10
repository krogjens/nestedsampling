function ns_test_evolve_plot(testlist)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% This routine plots in a single figure the results of ns_test_evolve_min
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

ncurves=length(testlist);
curves=zeros(length(testlist(1).res)+1,ncurves);
for i=1:ncurves
  curves(2:end,i)=testlist(i).res;
end
meanend=mean(curves(ceil(end/3):end,:));
for i=1:ncurves
  curves(:,i)=2*curves(:,i)/meanend(i)+i-1;
end

plot(transpose(0:(length(curves(:,1))-1)),curves)
hold on
plot(transpose(0:(length(curves(:,1))-1)),ones(length(curves(:,1)),1)*(2:(ncurves+1)),':')
hold off

