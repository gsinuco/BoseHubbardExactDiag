clear
clc


global D_H L N pi N_Floquet;

pi            = 4.0*atan(1.0);
L             = 5;
N             = 1;
D_H           = DimensionHilbert(L,N);
N_Floquet     = 1;                    % Number of Floquet channels = 2*N_Floquet + 1

N_phi         = 128;

%N_i           = zeros(L);
Eval          = zeros(N_phi,D_H);
Eig_Floquet   = zeros(N_phi,(2*N_Floquet+1)*D_H);
state_tag     = zeros(1,D_H);
stateK_tag    = zeros(1,D_H);
state         = zeros(D_H,L);
EigenEnergies = zeros(1,D_H);
H_Identity    = eye(D_H,D_H);
H_B           = sparse(D_H,D_H);
H_U           = sparse(D_H,D_H);
H_K           = sparse(D_H,D_H);
H_K           = sparse(D_H,D_H);
J_nn          = sparse(D_H,D_H);


[state, state_tag] = BuildingBasis();

t        = 1.0;        % Tunnelling factor
theta    = 0.0;%pi/L;%0.0;  % Gauge phase
lambda   = 0.0 ;
b        = 4.1;%(1.0 + sqrt(5.0))/2.0;

potencia = linspace(0,10,100);
U        = 1.0;
Sites    = linspace(1,L,L);
V        = zeros(1,L);

H_K = KineticEnergy(state, state_tag,t, theta);% + ...
      %InterparticleInteraction(state,state_tag,U);
phi = linspace(0,2*pi,N_phi);

for k_=1:size(phi,2)
    
    V      = lambda.*cos(b*2*pi*(Sites) + phi(k_));
    %U      = k_*2.0/size(phi,2);
    H      = H_K + InterparticleInteraction(state,state_tag,U) + LocalShifts(state,state_tag,V);
    H_full =  full(H);
    Eval(k_,1:D_H) = eig(H_full) ;
    
    
    %     omega     = 0.5; % Perturbation frequency
    %     V2        = 1.0; % second barrier height
    %
    %     H_F = zeros((2*N_Floquet+1)*D_H,(2*N_Floquet+1)*D_H);
    %
    %     for n = -N_Floquet:N_Floquet
    %         for m = -N_Floquet:N_Floquet
    %
    %             VF = 0.0;
    %             if n ~= m
    %                 for k=1:L
    %                     VF = VF + 1i*V2*exp(1i*(n-m)*k*2*pi/L)*(1 - exp(-1i*(n-m)*2*pi/L))/(2*pi*(n-m));
    %                 end
    %                 H_F((n+N_Floquet)*D_H + 1:(n+N_Floquet)*D_H + D_H,(m+N_Floquet)*D_H + 1:(m+N_Floquet)*D_H + D_H) = VF*H_Identity;
    %             else
    %                 H_F((n+N_Floquet)*D_H + 1:(n+N_Floquet)*D_H + D_H,(m+N_Floquet)*D_H + 1:(m+N_Floquet)*D_H + D_H) = H_full - n*omega*H_Identity;
    %             end
    %
    %         end
    %     end
    %
    %     Eig_Floquet(k_,1:(2*N_Floquet+1)*D_H) = eig(H_F);
    
end

%%
N_i    = zeros(L,D_H);
phi_   = 0.2*pi;
V      = lambda.*cos(b*2*pi*(Sites) + phi_);
H      = H_K +  LocalShifts(state,state_tag,V);
H_full =  full(H);
[V,D]  = eig(H_full) ;

%%% Local number of particles
N_i = 0.0;
for m_=1:D_H
    for l_=1:L
        n_l = diag(state(:,l_));       
        N_i(l_,m_) = V(:,m_)'*n_l*V(:,m_);
    end
end

%%% Current (Physica E 46, 119-132 (2012), Eq. (32))
J_nn   = Current(state, state_tag,t, theta, 1,2);
for m_=1:D_H
    J_full(m_) =  real(V(:,m_)'*full(J_nn)*V(:,m_));
end


%%
%%% Figures

figure(6)

subplot(2,2,1 )
plot(phi/pi,Eval)
hold on
%axis([0 2 -9.2 9.2])

for k_ = 1:1:D_H
    subplot(2,2,1)
    scatter(0.2,D(k_,k_))
    
    
    subplot(2,2,2)
    bar(N_i(:,k_))
    axis([0 L 0 0.7])
    hold off
    
    subplot(2,2,4)
    %plot(J_full);figure(gcf);
    scatter(k_,J_full(k_));
    axis([0 D_H     -0.1 0.1])
    hold on
    
    
    pause()
    
end

subplot(2,2,4)
plot(J_full);figure(gcf);
axis([0 D_H -0.1 0.1])
























%%
% figure(3)
% plot(phi,Eig_Floquet)
% axis([0 2*pi -14 14])
%%% tight-binding spectrum

%spy(H_full)
%figure(3)

%%E = TightBindingE(state,t,theta);