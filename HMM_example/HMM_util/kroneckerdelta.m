%implementation of kroneckerdelta without symbolic shitbox
%also works for doubles, so fuck matlab 
function z=kroneckerdelta(x,y)
    if x==y
        z=1;
    else
        z=0;
    end