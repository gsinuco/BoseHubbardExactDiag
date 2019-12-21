%% Kinetic Energy

function [H_Jx,H_Jy] = KineticEnergy2D(state,state_tag,alpha)


global D_H N L Lx Ly

H = 0.0;
state_k = 0;
 
k_Jx = 1;
k_Jy = 1;
for k_ = 1:Ly;
    for j_=1:Lx;
                
        %%% Boundary conditions - BEGIN
        if j_< Lx && k_<Ly
            jk   = (k_ - 1)*Lx    + j_;
            jkp1 = (k_ + 1 -1)*Lx + j_;
            jp1k = (k_ - 1)*Lx    + j_ + 1;
        end
        if j_==Lx && k_<Ly
            jk   = (k_ - 1)*Lx    + j_;
            jp1k = (k_ - 1)*Lx    + 1;
            jkp1 = (k_ + 1 -1)*Lx + 1;
        end
        if k_==Ly && j_<Lx
            jk   = (k_ - 1)*Lx + j_;
            jp1k = (k_ - 1)*Lx + j_ + 1;
            jkp1 = (1  - 1)*Lx + 1;
        end
        
        if k_==Ly && j_==Lx
            jk   = (k_ - 1)*Lx + j_;
            jp1k = (k_ - 1)*Lx + 1;
            jkp1 = (1  - 1)*Lx + 1;
        end
        
        %%% Boundary conditions - END
        
        for j = 1:D_H
            
            state_k = a_dagger(jk,a(jp1k,state(j,:)));
            
            if size(state_k,2)==L
              
                
                stateK_tag = stateTAG(state_k);
                I_Jx(k_Jx) = j;
              
                J_Jx(k_Jx) = find(state_tag==stateK_tag);
                s_Jx(k_Jx)  = sqrt((state(j,jk)+1)*state(j,jp1k));
                k_Jx = k_Jx + 1;
            end
            
            %             state_k = a_dagger(jk,a(jkp1,state(j,:)));
            %             if size(state_k,2)==L
            %                 stateK_tag = stateTAG(state_k);
            %                 I_(k) = j;
            %                 J_(k) = find(state_tag==stateK_tag);
            %                 s_(k)  = Jy*exp(-j_*theta*1i)*sqrt((state(j,jk)+1)*state(j,jkp1));
            %                 k = k + 1;
            %             end
            
            state_k = a_dagger(jk,a(jk,state(j,:))); % commented because
            %this term is diagonal
            if size(state_k,2)==L
                %stateK_tag = stateTAG(state_k);  % commented because this term is diagonal
                I_Jy(k_Jy) = j;
                J_Jy(k_Jy) = j;%find(state_tag==stateK_tag); %commented because this term is diagonal
                s_Jy(k_Jy)  = (cos(2*pi*j_*alpha + 2*pi*k_/Ly))*state(j,jk);
                k_Jy = k_Jy + 1;
            end

            
        end
        
        
    end
end

H_Jx = - sparse(I_Jx,J_Jx,s_Jx,D_H,D_H) - sparse(J_Jx,I_Jx,conj(s_Jx),D_H,D_H);
H_Jy = - sparse(I_Jy,J_Jy,s_Jy,D_H,D_H) - sparse(J_Jy,I_Jy,conj(s_Jy),D_H,D_H);
%H =  H - sparse(I_,J_,s_Jx,D_H,D_H) - sparse(J_,I_,conj(s_Jx),D_H,D_H) - sparse(I_,J_,s_Jy,D_H,D_H) - sparse(J_,I_,conj(s_Jy),D_H,D_H);

end