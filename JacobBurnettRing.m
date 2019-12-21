clear
clc


global L N pi;

pi  = 4.0*atan(1.0);
L   = 4;
N   = 4;



D_H = DimensionHilbert(L,N);

Eval = zeros(100,D_H);
state_tag = zeros(1,D_H);
stateK_tag= zeros(1,D_H);

state(1,:) = zeros(1,L);
state(1,1) = N;

state_tag(1) = state(1,1);

state(2,:) = zeros(1,L);
state(2,1) = N-1;
state(2,2) = 1;
state_tag(2) = state(2,1) + state(2,2)*(N+1);

%% Building the basis

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

%% Hamiltonian Parameters

EigenEnergies = zeros(1,D_H);
H_UandB       = sparse(D_H,D_H);
H             = sparse(D_H,D_H);
H_aux         = sparse(D_H,D_H);

%% Barrier
V     = zeros(1,L);
V(2) = 0.0;
V(3) = 0.0;
V(1) = 0.0;

I = linspace(1,D_H,D_H);
J = linspace(1,D_H,D_H);

s = 0;
for j=1:D_H
    s(j) = 0;
    for k=1:L
        s(j) = s(j) + V(k)*state(j,k);
    end
end

H_UandB =  sparse(I,J,s,D_H,D_H);

%% on-site interaction

potencia = linspace(-4,4,100);
U = 0.50;
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
H_UandB = H_UandB + sparse(I,J,s,D_H,D_H);

%% Kinetic Energy

t = 1.0;
t_l = 1.0;
t_r = 1.0;
theta = linspace(0,2*pi/L,100);
E = zeros(size(theta,2),D_H);

for m=1:100;
    H = 0.0;
    %s = 0.0;
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
                
                I_(k) = j;
                J_(k) = find(state_tag==stateK_tag(j));
                if i<L
                    s_(k) = sqrt((state(j,i) + 1)*state(j,i+1));
                else
                    s_(k) = sqrt((state(j,i) + 1)*state(j,1));
                end
                k = k + 1;
                
            end
        end
        M = mod(i,2);
        if M==0
            t = t_l;
        else
            t= t_r;
        end
        H =  H - t*exp(-theta(m)*1i)*sparse(I_,J_,s_,D_H,D_H) - t*exp(theta(m)*1i)*sparse(J_,I_,s_,D_H,D_H);        
    
    end
    
 H= H + H_UandB;
 

%% tight-binding spectrum


for j=1:D_H
    for q=0:L-1
       E(m,j) =  E(m,j) - 2*t*cos(2*q*pi/L -theta(m))*state(j,q+1);
    end
end

%%

%[Evec,Eval] = eigs(H,3,'sa');
H_full =  full(H);
Eval(m,1:D_H) = eig(H_full) ;
%[Evec,E_] = eig(H_full) ;

end

E = sort(E,2);


%%
figure(2)
spy(H_full)
figure(3)
plot(theta/(2*pi/L),Eval);
figure(4)

plot(theta/(2*pi/L),E);

% 



