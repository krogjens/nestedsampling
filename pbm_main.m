function out = pbm_main(obs,mu,sigma,out_wish)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Contributors to the programming: Michael Lomholt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if out_wish=='l'
  var=sigma^2;
  stepdevs=obs-mu(1,:);
  out=-log(2*pi*var)*numel(stepdevs)/2-sum(sum(stepdevs.^2))/(2*var);
elseif out_wish=='u'
  stepdevs=obs-mu(1,:);
  out=1-erfc(stepdevs/(sqrt(2)*sigma))/2;
elseif out_wish=='x'
  stepdevs = sqrt(2)*sigma*erfcinv(2-2*obs);
  out = stepdevs+mu(1,:);
end
end
