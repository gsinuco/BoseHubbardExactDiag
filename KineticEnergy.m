%% Kinetic Energy

function [H] = KineticEnergy(state,state_tag,t,theta)


global D_H N L

H = 0.0;
state_k = 0;

for i = 1:L-1;
    
    k = 1;
    
    for j = 1:D_H
        
        if i<L
            state_k = a_dagger(i,a(i+1,state(j,:)));
        else
            
            state_k = a_dagger(i,a(1,state(j,:))); % periodic boundary conditions
        end
        %state_k;
        if size(state_k,2)==L
            
            stateK_tag = stateTAG(state_k);
            
             I_(k) = j;
%             stateK_tag;
%             find(state_tag==stateK_tag)
            J_(k) = find(state_tag==stateK_tag);
            if i<L
                s_(k) =  sqrt((state(j,i) + 1)*state(j,i+1));
            else
                s_(k) =  sqrt((state(j,i) + 1)*state(j,1));
            end
            k = k + 1;
            
        end
    end
    
    H =  H - t*exp(-theta*1i)*sparse(I_,J_,s_,D_H,D_H) - t*exp(theta*1i)*sparse(J_,I_,s_,D_H,D_H);
    
end