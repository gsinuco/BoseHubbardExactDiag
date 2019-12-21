clear
clc

global D_H L N pi N_Floquet Lx Ly;

pi            = 4.0*atan(1.0);
Lx            = 64;
Ly            = 64;
L             = Lx*Ly;
N             = 1;
D_H           = DimensionHilbert(L,N);
N_Floquet     = 1;

%Eval          = zeros(N_phi,D_H);
%Eig_Floquet   = zeros(N_phi,(2*N_Floquet+1)*D_H);
state_tag     = zeros(1,D_H);
stateK_tag    = zeros(1,D_H);
state         = zeros(D_H,L);
EigenEnergies = zeros(1,D_H);
H_Identity    = eye(D_H,D_H);
%H_B           = sparse(D_H,D_H);
%H_U           = sparse(D_H,D_H);
%H_K           = sparse(D_H,D_H);
%J_nn          = sparse(D_H,D_H);
H_kyIdentity  = eye(Lx,Lx);
H_JxS         = sparse(D_H,D_H);
H_JyS         = sparse(D_H,D_H); 

[state, state_tag] = BuildingBasis();

Jx     = pi;
Jy     = pi;
alpha  = 1.0/3.0;

%%% A)
tau = 1.0/3.0;
T   = 5.0/6.0;


omega = pi/T;

%%% B)
Jx_A = 1.5;
Jy_A = 1.5;
phi  = pi/2.0;


U      = 1.0;
[H_JxS,H_JyS] = KineticEnergy2D(state, state_tag,alpha);
H      = Jx*H_JxS + Jy*H_JyS;
H_full = full(H);
H_Jx   = full(H_JxS);
H_Jy   = full(H_JyS);
%[V,D]  = eig(H_full) ;

%%% Spectrum of the static system
for ky=-Ly:-1
    H_ky = H_full(1+(ky+Ly)*Lx:(ky+Ly+1)*Lx,1+(ky+Ly)*Lx:(ky+Ly+1)*Lx);
    [V_ky,D_ky]  = eig(H_ky) ;
    E(:,ky+Ly+1) = diag(D_ky);
    
    
    H_F    = zeros((2*N_Floquet+1)*Lx,(2*N_Floquet+1)*Lx);
    Hky_Jx = H_Jx(1+(ky+Ly)*Lx:(ky+Ly+1)*Lx,1+(ky+Ly)*Lx:(ky+Ly+1)*Lx);
    Hky_Jy = H_Jy(1+(ky+Ly)*Lx:(ky+Ly+1)*Lx,1+(ky+Ly)*Lx:(ky+Ly+1)*Lx);
    for n = -N_Floquet:N_Floquet
        for m = -N_Floquet:N_Floquet
            
            
            if n ~= m
                %%% A)
                VF_Jx = Jx*1i*(exp(-1i*(m-n)*omega*tau) - 1)/(2*pi*(m-n));
                VF_Jy = Jy*1i*(exp(-1i*(m-n)*omega*(T-tau)) - exp(-1i*(m-n)*omega*tau))/(2*pi*(m-n));
                %%% B)
                %if(abs(m-n)==1)
                %    VF_Jx = 0.5*Jx_A;
                %    VF_Jy = 0.5*Jy_A*exp(1i*(m-n)*phi);
                %else
                %    VF_Jx = 0.0;
                %    VF_Jy = 0.0;
                %end
                H_F((n+N_Floquet)*Lx + 1:(n+N_Floquet)*Lx + Lx,(m+N_Floquet)*Lx + 1:(m+N_Floquet)*Lx + Lx) = VF_Jx*Hky_Jx + VF_Jy*Hky_Jy;
            else
                %%% A)
                H_F((n+N_Floquet)*Lx + 1:(n+N_Floquet)*Lx + Lx,(m+N_Floquet)*Lx + 1:(m+N_Floquet)*Lx + Lx) = Jx*(tau)*Hky_Jx + Jy*(T-tau)*Hky_Jy - n*omega*H_kyIdentity;
                %%% B)
                %H_F((n+N_Floquet)*Lx + 1:(n+N_Floquet)*Lx + Lx,(m+N_Floquet)*Lx + 1:(m+N_Floquet)*Lx + Lx) = Jx*Hky_Jx + Jy*Hky_Jy - n*omega*H_kyIdentity;
            end
            
        end
    end
    
    Eig_Floquet(:,ky+Ly+1) = eig(H_F);
end


%%
figure(1)
for ky=-Ly:-1
    hold on
    plot(E(ky+Ly+1,:))
end


%%
figure(3)
for ky=1:size(Eig_Floquet,1)
    hold on
    plot(Eig_Floquet(ky,:))  
end