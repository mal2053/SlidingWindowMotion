function [V] = Matrix2Vector(M,Q)    
    V = M(Q==1);
end