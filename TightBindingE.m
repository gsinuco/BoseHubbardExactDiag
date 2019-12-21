%% tight-binding spectrum
function [E] = TightBindingE(state,t,theta)

global D_H N L pi

E = zeros(size(theta,2),D_H);
for j=1:D_H
    for q=0:L-1
       E(j) =  E(j) - 2*t*cos(2*q*pi/L -theta)*state(j,q+1);
    end
end
E = sort(E,2);