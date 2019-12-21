
clear
clc


global D_H L N pi;

pi  = 4.0*atan(1.0);
L   = 20;
N   = 1;
D_H = DimensionHilbert(L,N);
Eval          = zeros(100,D_H);
state_tag     = zeros(1,D_H);
stateK_tag    = zeros(1,D_H);
state         = zeros(D_H,L);
EigenEnergies = zeros(1,D_H);
H_Identity    = eye(D_H,D_H);
H_B           = sparse(D_H,D_H);
H_U           = sparse(D_H,D_H);
H_K           = sparse(D_H,D_H);
[state, state_tag] = BuildingBasis();



t        = 1.0;
theta    = 0.0;
potencia = linspace(0,10,100);
U        = linspace(0,10,100);
V        = zeros(1,L);
V(1)     = 0.0;
V(2)     = 0.0;
V(3)     = 0.0;
V(4)     = 0.0;
V(5)     = 0.0;



H_K = KineticEnergy(state, state_tag,t, theta);
H_K = H_K +  LocalShifts(state,state_tag,V);
for m=1:size(potencia,2)
    
    H = H_K + InterparticleInteraction(state,state_tag,U(m));
          
    %[Evec,Eval] = eigs(H,3,'sa');
    H_full =  full(H);
    Eval(m,1:D_H) = eig(H_full) ;
    %[Evec,E_] = eig(H_full) ;
    
end

%% tight-binding spectrum
E = TightBindingE(state,t,theta);


%%
figure(2)
spy(H_full)
figure(3)
plot(U,Eval)
axis([0 10 -14 14])



% %% FLOQUET
% 
% N_Floquet = 5.0; % Number of Floquet channels = 2*N_Floquet + 1
% omega     = 0.5; % Perturbation frequency
% V2        = 1.0; % second barrier height
% 
% H_F = zeros((2*N_Floquet+1)*D_H,(2*N_Floquet+1)*D_H);
% 
% for n = -N_Floquet:N_Floquet
% for m = n:N_Floquet
%     
%     
%     VF = 0.0;
%     for k=1:L
%         VF = VF + 1i*V2*exp(1i*(n-m)*k*2*pi/L)*(1 - exp(-1i*(n-m)*2*pi/L));
%     end
%     
%     H_F((n+N_Floquet)*D_H + 1:(n+N_Floquet)*D_H + 1+D_H,m+N_Floquet + 1:m+N_Floquet + 1+D_H) = VF*H_Identity + (H_full - n*omega)*kroneckerDelta(n,m);
% end
% end

  
  
