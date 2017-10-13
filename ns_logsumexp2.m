function z = logsumexp2(x, y)
% z = log(exp(x)+exp(y))

if x == -Inf && y == -Inf
   z = -Inf;
else
   if x>y
     z = x+log(1+exp(y-x));
   else
     z = y+log(1+exp(x-y));
   end
end

%log(A)+log(1+exp(log(B)-log(A)))
%log(A)+log(1+exp(log(B/A))
%log(A)+log(1+B/A))
%log(A(1+B/A))
%log(A+B) QED
