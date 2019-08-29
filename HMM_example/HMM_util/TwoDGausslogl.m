%2-D gaussian likelihood function rewritten to log(f(x)) 
%and further rewrite to pull log(2pi) out of computation loop to save a bit
%of time
%mu contains mean for x and y dimension respectively
%sigma is the standard deviation


function b = TwoDGausslogl(x,y,sigma, mu_x,mu_y, log2pi)
     b=((x-mu_x).*(mu_x-x)+(y-mu_y).*(mu_y-y))./(sigma^2*2)-2*log(sigma)-log2pi;
end
