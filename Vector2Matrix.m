function [M] = Vector2Matrix(V,Q)    
    M = zeros(size(Q));
    M(Q == 1) = V;
    M = M + M';
end