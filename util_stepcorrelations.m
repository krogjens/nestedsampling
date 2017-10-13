function [vals]=util_stepcorrelations(obs,nsteps)
% Estimates the correlations between single steps seperated by the values in nsteps

vals=zeros(1,length(nsteps));
len=length(obs)-1;
steps=obs(2:(len+1))-obs(1:len);
for i=1:length(nsteps)
  for j=1:(len-nsteps(i))
    vals(i)=vals(i)+steps(j)*steps(j+nsteps(i));
  end
  vals(i)=vals(i)/(len-nsteps(i));
end

