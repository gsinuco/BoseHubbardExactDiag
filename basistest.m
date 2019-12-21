clear
clc

L = 5;
N = 3;

state(1,:)   = zeros(1,L+2);
state(1,1) = N;
state(1,2) = 0;
state_alt(1,:)   = zeros(1,L+2);
state_alt(1,1) = N;
state_alt(1,2) = 0;

k = 0;
j = 1;
count = 0
offset = 0
for l = 2:13;%16
     
    
    state_alt(l)
    while k<L-1
        N_k = 0;
        for i=1:k+1
            N_k = N_k + state(j,i);
        end
        if k> 0 
            state(j+1,1:k)   = state(j,1:k);
        end
        state(j+1,k+1)   = state(j,k+1) - 1;
        state(j+1,k+2)   = N - (N_k - 1);
        state(j+1,k+3:L) = 0;
        state(j+1,L+1)   = l;
        state(j+1,L+2) = k;                             
        
        k = k+1;        
        j = j+1;
        
    end
    [C,minpos] = min(state(j,:));
    
    if(state(j,L) == N-1) 
        offset = offset + 1
    end
    l
    offset
    if(count > 0) 
     k = mod(l + 2 + offset -1,L-1) ;
    else
     k = mod(l + 2,L-1) ;
    end
    
    
    %while
    
    
%     if l == 2 
%         k = 0;
%     end
%     if l == 3
%         k = 1;
%     end
%         
%     if l == 4
%         k = 2;
%     end
%     
%     if l == 5
%         k = 3;
%     end
%     
%     if l == 6
%         k = 0;
%     end
%     
%     if l == 7 
%         k = 1;
%     end
%     if l == 8
%         k = 2;
%     end
%         
%     if l == 9
%         k = 3;
%     end
%     
    if l == 10
        k=1;
    end
     
    if l == 11
        k=2;
    end
    
     
    if l == 12
        k=3;
    end
    
     
    if l == 13
        k=2;
    end
    
    k
%     if l == 14
%         k=3;
%     end
%     if l == 15
%         k=3;
%     end
end
    
    
 
%end
