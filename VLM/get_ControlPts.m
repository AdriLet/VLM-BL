function [Sail] = get_ControlPts (Sail)

X_V=Sail.X_V;
Y_V=Sail.Y_V;
Z_V=Sail.Z_V;

M=length(X_V(:,1));
N=length(X_V(1,:));

for j=1:N-1   
    for i=1:M-1
        X_C(i,j)= (X_V(i,j)+X_V(i+1,j)+X_V(i+1,j+1)+X_V(i,j+1))/4;
        Y_C(i,j)= (Y_V(i,j)+Y_V(i+1,j)+Y_V(i+1,j+1)+Y_V(i,j+1))/4;
        Z_C(i,j)= (Z_V(i,j)+Z_V(i+1,j)+Z_V(i+1,j+1)+Z_V(i,j+1))/4;         
    end      
end

Sail.X_C=X_C;
Sail.Y_C=Y_C;
Sail.Z_C=Z_C;