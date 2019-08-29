function z = LOGPLUS(x,y) %does log(exp(x)+exp(y))
    if x>y
       z= x+log(1+exp(y-x));
    else
       z= y+log(1+exp(x-y));
    end
end