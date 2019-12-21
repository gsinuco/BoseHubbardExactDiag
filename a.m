function [state_] = a(i,state)

global L
state_ = zeros(1,L);
if(state(i)>0)
    if(i>=2) 
        state_(1:i-1) = state(1:i-1);
    end
    state_(i)     = state(i)-1;
    if i<L
        state_(i+1:L)     = state(i+1:L);
    end
else
    state_ = 0;
end



end