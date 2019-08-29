%%hmm data generator 2D diffusion n states
%delta x is data 
function [dx,dy,stateseq] = HMMgenerator(p,A,B,u,steps) %inverse pdf also or some shit later
%transition matrix between states 

    %inital state
    state=sum(rand >= cumsum(p))+1;
    dx=zeros(1,steps);
    dy=zeros(1,steps);
    stateseq=zeros(1,steps);
    for h = 1:steps 
        %emission
        dx(h)=normrnd(u(state,1), B(state),1,1);
        dy(h)=normrnd(u(state,2), B(state),1,1);
        %state change
        stateseq(h)=state;
        state=sum(rand >= cumsum(A(state,:)))+1;
    end
end