function z = logsumexp2(x, y)
if x>y
  z = x+log(1+exp(y-x));
else
  z = y+log(1+exp(x-y));
end
