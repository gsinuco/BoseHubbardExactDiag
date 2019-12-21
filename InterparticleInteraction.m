%% Interparticle interaction

function [H_U] = InterparticleInteraction(state,state_tag,U)

global D_H L

I = linspace(1,D_H,D_H);
J = linspace(1,D_H,D_H);
s = 0.0 ;
for j=1:D_H
    s(j) = 0;
    for k=1:L
        s(j) = s(j) + state(j,k)*(state(j,k)-1);
    end
end
s = 0.5*U*s;

H_U = sparse(I,J,s,D_H,D_H);
