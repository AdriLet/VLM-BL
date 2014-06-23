function [Sail] = get_normal (Sail)

M = length(Sail.X_V(:,1));
N = length(Sail.X_V(1,:));

for i=1:M-1   
    for j=1:N-1       
        x1=Sail.X_V(i,j);
        y1=Sail.Y_V(i,j);
        z1=Sail.Z_V(i,j);
        
        x2=Sail.X_V(i+1,j);
        y2=Sail.Y_V(i+1,j);
        z2=Sail.Z_V(i+1,j);
        
        x3=Sail.X_V(i,j+1);
        y3=Sail.Y_V(i,j+1);
        z3=Sail.Z_V(i,j+1);
        
        x4=Sail.X_V(i+1,j+1);
        y4=Sail.Y_V(i+1,j+1);
        z4=Sail.Z_V(i+1,j+1);

        %Definition of the Vectors A and B

        Ax=x4-x1;
        Ay=y4-y1;
        Az=z4-z1;

        Bx=x3-x2;
        By=y3-y2;
        Bz=z3-z2;

        % The vector unit normal

        nx (i,j)= (Ay*Bz-By*Az )/ ((Ay*Bz-By*Az)^2 + (-Ax*Bz+Bx*Az)^2 + (Ax*By-Bx*Ay)^2 )^0.5;
        ny (i,j)= (-Ax*Bz+Bx*Az) / ((Ay*Bz-By*Az)^2 + (-Ax*Bz+Bx*Az)^2 + (Ax*By-Bx*Ay)^2 )^0.5;
        nz (i,j)= (Ax*By-Bx*Ay) / ((Ay*Bz-By*Az)^2 + (-Ax*Bz+Bx*Az)^2 + (Ax*By-Bx*Ay)^2 )^0.5;
        
    end   
end

Sail.nx=nx;
Sail.ny=ny;
Sail.nz=nz;