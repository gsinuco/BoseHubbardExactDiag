%% Current

function [J] = Current(state,state_tag,t,theta,i,j)


global D_H L

J = 0.0;
state_k = 0;
    
k = 1;
for m = 1:D_H
    
    state_k = a_dagger(i,a(j,state(m,:)));    

    if size(state_k,2)==L
        
        stateK_tag = stateTAG(state_k);
        I_(k) = m;
        J_(k) = find(state_tag==stateK_tag);
        s_(k) = sqrt(state(m,i)+1)*sqrt(state(m,j));
        
        k = k + 1;
        
    end
end

J =  -1i*(t*exp(-theta*1i)*sparse(I_,J_,s_,D_H,D_H) - t*exp(theta*1i)*sparse(J_,I_,s_,D_H,D_H));

%end