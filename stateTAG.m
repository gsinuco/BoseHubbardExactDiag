function [ST]=stateTAG(state)

global L N

ST = 0;
% for l=1:L
%     ST = ST + state(l)*(N+1)^(l-1);
% end
% 
for l=1:L
    p_l = 100.0*l+3.0;
    ST = ST + sqrt(p_l)*state(l);
end
        
end