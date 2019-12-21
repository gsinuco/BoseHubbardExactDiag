function [D_H] = DimensionHilbert(L,N)

%L = 5
%N = 3

n_ = L+N-1;
k_ = N;

Num = 1.0;
for i = 1:k_
    Num = Num*n_;
    n_ = n_-1;
end
Den = 1.0;
for i = 1:k_-1
    Den = Den*k_;
    k_ = k_-1;
end

D_H = Num/Den;

end