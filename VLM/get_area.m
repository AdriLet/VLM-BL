function [Sail] = get_area (Sail)

M =length(Sail.X_V(:,1));
N =length(Sail.X_V(1,:));

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
        
        %Definition of the Vectors C and D
        Cx=x3-x1;
        Cy=y3-y1;
        Cz=z3-z1;

        Dx=x2-x1;
        Dy=y2-y1;
        Dz=z2-z1;

        % The dot product
        nx = (Cy*Dz-Dy*Cz );
        ny = (-Cx*Dz+Dx*Cz);
        nz = (Cx*Dy-Dx*Cy);

        % The Area (=mangitude)
        Sail.area (i,j)= (nx^2+ny^2+nz^2)^0.5;        
    end    
end