clear
clc

global D_H L N pi N_Floquet Lx Ly;

pi            = 4.0*atan(1.0);
Lx            = 32;
Ly            = 32;
L             = Lx*Ly;
N             = 2;
D_H           = DimensionHilbert(L,N);
N_Floquet     = 1;

%Eval          = zeros(N_phi,D_H);
%Eig_Floquet   = zeros(N_phi,(2*N_Floquet+1)*D_H);
state_tag     = zeros(1,D_H);
% stateK_tag    = zeros(1,D_H);
state         = zeros(D_H,L);
% EigenEnergies = zeros(1,D_H);
% H_Identity    = eye(D_H,D_H);
% H_B           = sparse(D_H,D_H);
% H_U           = sparse(D_H,D_H);
% H_K           = sparse(D_H,D_H);
% J_nn          = sparse(D_H,D_H);
% H_Jx          = sparse(D_H,D_H);
% H_Jy          = sparse(D_H,D_H); 


[state, state_tag] = BuildingBasis();
% 
% Jx     = 1.0;
% Jy     = 1.0;
% alpha  = 1.0/3.0;
% 
% U           = 1.0;
% [H_Jx,H_Jy] = KineticEnergy2D(state, state_tag,alpha);
% H      = Jx*H_Jx + Jy*H_Jy;
% H_full =  full(H);
% %[V,D]  = eig(H_full) ;
% 
% for ky=-Ly:-1     
%     H_ky = H_full(1+(ky+Ly)*Lx:(ky+Ly+1)*Lx,1+(ky+Ly)*Lx:(ky+Ly+1)*Lx);
%     [V_ky,D_ky]  = eig(H_ky) ;
%     E(:,ky+Ly+1) = diag(D_ky);
% end
% 
% %%
% figure(2)
% for ky=-Ly:-1
%     hold on
%     plot(E(ky+Ly+1,:))  
% end
% 
% %%% Local number of particles
% N_i = 0.0;
% for m_=1:D_H
%     for l_=1:L
%         n_l = diag(state(:,l_));       
%         N_i(l_,m_) = V(:,m_)'*n_l*V(:,m_);
%     end
% end
% 
% 
% figure(6)
% 
% for k_ = 1:1:D_H
%     subplot(2,2,1)
%     scatter(0.2,D(k_,k_))    
%     subplot(2,2,2)
%     bar(N_i(:,k_))
%     axis([0 L 0 0.7])
%     hold off    
%     subplot(2,2,4)
%     scatter(k_,D(k_,k_));
%     hold on        
% %    pause()    
% end
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% 
% 
% %%% Current (Physica E 46, 119-132 (2012), Eq. (32))
% % J_nn   = Current(state, state_tag,t, theta, 1,2);
% % for m_=1:D_H
% %     J_full(m_) =  real(V(:,m_)'*full(J_nn)*V(:,m_));
% % end
% 
% 
% % + ...
%       %InterparticleInteraction(state,state_tag,U);
% %phi = linspace(0,2*pi,N_phi);
% % 
% % for k_=1:1;%size(phi,2)
% %     
% %     %U      = k_*2.0/size(phi,2);
% %     H      = H_K + 0.0*InterparticleInteraction(state,state_tag,U) + 0.0*LocalShifts(state,state_tag,V);
% %     H_full =  full(H);
% %     Eval(k_,1:D_H) = eig(H_full) ;
% %     
% %     
% %     %     omega     = 0.5; % Perturbation frequency
% %     %     V2        = 1.0; % second barrier height
% %     %
% %     %     H_F = zeros((2*N_Floquet+1)*D_H,(2*N_Floquet+1)*D_H);
% %     %
% %     %     for n = -N_Floquet:N_Floquet
% %     %         for m = -N_Floquet:N_Floquet
% %     %
% %     %             VF = 0.0;
% %     %             if n ~= m
% %     %                 for k=1:L
% %     %                     VF = VF + 1i*V2*exp(1i*(n-m)*k*2*pi/L)*(1 - exp(-1i*(n-m)*2*pi/L))/(2*pi*(n-m));
% %     %                 end
% %     %                 H_F((n+N_Floquet)*D_H + 1:(n+N_Floquet)*D_H + D_H,(m+N_Floquet)*D_H + 1:(m+N_Floquet)*D_H + D_H) = VF*H_Identity;
% %     %             else
% %     %                 H_F((n+N_Floquet)*D_H + 1:(n+N_Floquet)*D_H + D_H,(m+N_Floquet)*D_H + 1:(m+N_Floquet)*D_H + D_H) = H_full - n*omega*H_Identity;
% %     %             end
% %     %
% %     %         end
% %     %     end
% %     %
% %     %     Eig_Floquet(k_,1:(2*N_Floquet+1)*D_H) = eig(H_F);
% %     
% % end
