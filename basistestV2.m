clear
clc


global L N pi;



pi  = 4.0*atan(1.0);
L   = 5;
N   = 3;

t = 1.0;
U = 1.0*t;
theta = pi/L;

D_H = DimensionHilbert(L,N);

H = sparse(D_H,D_H);

state_tag = zeros(1,D_H);
stateK_tag= zeros(1,D_H);

state(1,:) = zeros(1,L);
state(1,1) = N;

state_tag(1) = state(1,1);

state(2,:) = zeros(1,L);
state(2,1) = N-1;
state(2,2) = 1;
state_tag(2) = state(2,1) + state(2,2)*(N+1);



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
        k_bar = I_max - (nonzerosat(1,size(nonzerosat,2)) - nonzerosat(1,size(nonzerosat,2)));
        state(j+1,1:k_bar) = state(j,1:k_bar);
        state(j+1,I_max)   = state(j,I_max) - 1;
        state(j+1,I_max+1) = state(j,I_max+1) + 1;
    end
    
    state_tag(j+1) = stateTAG(state(j+1,:));
end


for i = 1:L;
    
    k = 1;
    
    for j = 1:D_H
        
        if i<L            
            state_k = a_dagger(i,a(i+1,state(j,:)));
        else
            state_k = a_dagger(i,a(1,state(j,:)));
        end
        
        if size(state_k,2)==L
           
            stateK_tag(j) = stateTAG(state_k);
            
            I(k) = j;
            J(k) = find(state_tag==stateK_tag(j));
            if i<L 
                 s(k) = sqrt((state(j,i) + 1)*state(j,i+1));
            else
                 s(k) = sqrt((state(j,i) + 1)*state(j,1));
            end
            k = k + 1;
            
        end
    end
    
    H = H  - t*exp(theta*1i)*sparse(I,J,s,D_H,D_H) - t*exp(-theta*1i)*sparse(J,I,s,D_H,D_H);
    
end

I = linspace(1,D_H,D_H);
J = linspace(1,D_H,D_H);

s = 0;
for j=1:D_H
    s(j) = 0;
    for k=1:L
        s(j) = s(j) + state(j,k)*(state(j,k)-1);
    end
end
s = 0.5*U*s;
H =  H + sparse(I,J,s,D_H,D_H);

figure(2)
spy(H)

%[Evec,Eval] = eigs(H,3,'sa');
[Evec,Eval] = eigs(H);
