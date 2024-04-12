function [C, Theta] = DCC(Dat)
% function [Ct, Ht] = DCC(dat)

[T,p] = size(Dat);
C = ones(p,p,T);
Theta = zeros(p,p,8);

for j=1:(p-1),    
    for k=(j+1):p,
  
        [ Ctmat , ~, param] = DCCsimple([Dat(:,j) Dat(:,k)]); 
        C(j,k,:) = Ctmat(1,2,:); 
        C(k,j,:) = Ctmat(1,2,:); 

        Theta(j,k,:) = param;
    end
end


