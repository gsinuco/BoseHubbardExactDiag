function [state, state_tag] = BuildingBasis()

global D_H N L
% state(j,:) : ocupation number basis. 
state(1,:)   = zeros(1,L);
state(1,1)   = N;

%state_tag(1) = state(1,1);
% ST = 0.0;
% for l=1:L
%     p_l = 200.0*l+3.0;
%     ST = ST + sqrt(p_l)*state(1,l);
% end
ST = stateTAG(state(1,:));
state_tag(1) = ST;

state(2,:)   = zeros(1,L);
state(2,1)   = N-1;
state(2,2)   = 1;
%state_tag(2) = state(2,1) + state(2,2)*(N+1);
%ST = 0.0;
% for l=1:L
%     p_l = 200.0*l+3.0;
%     ST = ST + sqrt(p_l)*state(2,l);
% end
ST = stateTAG(state(2,:));
state_tag(2) = ST;

k_bar        = 1; 

for j = 2:D_H-1
    
    zerosat    = find(state(j,:)==0);
    nonzerosat = find(state(j,:));
    I_max      = find(state(j,:),1,'last');
    
    if size(nonzerosat,2)>1
        k_bar = I_max - (nonzerosat(1,size(nonzerosat,2)) - nonzerosat(1,size(nonzerosat,2)-1));
    end
    
    N_k =  0;
    for i=1:k_bar
        N_k =  N_k + state(j,i);
    end
    
    if (state(j,I_max-1)>0 && I_max+1<=L)
        state(j+1,1:I_max-1) = state(j,1:I_max-1);
        state(j+1,I_max)     = state(j,I_max) - 1;
        state(j+1,I_max+1)   = state(j,I_max+1) + 1;
    end
    
    if (I_max == L && state(j,I_max-1)==0)
        k_bar = I_max - (nonzerosat(1,size(nonzerosat,2)) - nonzerosat(1,size(nonzerosat,2)-1));
        state(j+1,1:k_bar-1) = state(j,1:k_bar-1);
        state(j+1,k_bar)     = state(j,k_bar) -1;
        state(j+1,k_bar + 1) = N - (N_k -1);
    end
    
    if (I_max == L && state(j,I_max-1)>0);
        state(j+1,1:I_max-2) = state(j,1:I_max-2);
        state(j+1,I_max-1)   = state(j,I_max-1) -1;
        state(j+1,I_max)     = state(j,I_max) +1;
    end
    
    if (I_max < L && state(j,I_max-1) == 0)
        k_bar = I_max ;
        state(j+1,1:k_bar) = state(j,1:k_bar);
        state(j+1,I_max)   = state(j,I_max) - 1;
        state(j+1,I_max+1) = state(j,I_max+1) + 1;
    end
    
    state_tag(j+1) = stateTAG(state(j+1,:));
end