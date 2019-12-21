    %% Barrier
function [H_B] = LocalShifts(state,state_tag,V)

global D_H N L

I = linspace(1,D_H,D_H);
J = linspace(1,D_H,D_H);

s = 0;
for j=1:D_H
    s(j) = 0;
    for k=1:L
        s(j) = s(j) + V(k)*state(j,k);
    end
end

H_B =  sparse(I,J,s,D_H,D_H);
    