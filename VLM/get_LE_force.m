function [Sail] = get_LE_force (Sail,rho,Le)

if Le==0
    
    FLeX=0;
    FLeY=0;
    FLeZ=0;
    MLeX=0;
    MLeY=0;
    MLeZ=0;
    
end
   
if Le==1

    X_V=Sail.X_V;
    Y_V=Sail.Y_V;
    Z_V=Sail.Z_V;

    X_C=Sail.X_C;
    Y_C=Sail.Y_C;
    Z_C=Sail.Z_C;

    nx= Sail.nx;
    ny= Sail.ny;
    nz= Sail.nz;

    GAMMA=Sail.GAMMA;

    N=length(Sail.Y_V(1,:));


    % Leading edge force, Guermond

    for j=1:N-1
       
        % Computation of the shortest distance from the CP to the LE (S)
        Ux=X_C(1,j)-X_V(1,j);
        Uy=Y_C(1,j)-Y_V(1,j);
        Uz=Z_C(1,j)-Z_V(1,j);
   
        Vx=X_V(1,j+1)-X_V(1,j);
        Vy=Y_V(1,j+1)-Y_V(1,j);
        Vz=Z_V(1,j+1)-Z_V(1,j);
   
        [Proj_UontoV,Proj_UontoV_x,Proj_UontoV_y,Proj_UontoV_z] = get_projection (Ux,Uy,Uz,Vx,Vy,Vz);
   
        S(1,j)= ((Ux-Proj_UontoV_x)^2 + (Uy-Proj_UontoV_y)^2 + (Uz-Proj_UontoV_z)^2 )^0.5;
   
   
        % Computation of the Force (it's a scalar)
   
        Lef_coef(1,j) = GAMMA(1,j) / (S(1,j)^0.5);
        LE_force_perunit(1,j) = pi*rho*(Lef_coef(1,j)^2) / 16;
   
   
        li(1,j)=[(X_V(1,j+1) - X_V(1,j))^2 + (Y_V(1,j+1) - Y_V(1,j))^2 + (Z_V(1,j+1) - Z_V(1,j))^2] ^0.5;
        LE_force(1,j) = LE_force_perunit (1,j)* li(1,j);
   
        % Computation of the Unit vector tangent of the Leading Edge
        [tX(1,j),tY(1,j),tZ(1,j)] = get_UnitVect (X_V(1,j),Y_V(1,j),Z_V(1,j),X_V(1,j+1),Y_V(1,j+1),Z_V(1,j+1));
    
        fiX(1,j) = ny(1,j)*tZ(1,j) - tY(1,j)*nz(1,j);
        fiY(1,j) = - (nx(1,j)*tZ(1,j) - tX(1,j)*nz(1,j));
        fiZ(1,j) = nx(1,j)*tY(1,j) - (tX(1,j)*ny(1,j));
   
        FiX(1,j) = fiX(1,j) / (fiX(1,j)^2+fiY(1,j)^2+fiZ(1,j)^2)^0.5;
        FiY(1,j) = fiY(1,j) / (fiX(1,j)^2+fiY(1,j)^2+fiZ(1,j)^2)^0.5;
        FiZ(1,j) = fiZ(1,j) / (fiX(1,j)^2+fiY(1,j)^2+fiZ(1,j)^2)^0.5;

   
        % Projection of the force onto the Unit vector perpendicular to the LE
        FLeX (1,j) = LE_force(1,j) * FiX(1,j);
        FLeY (1,j) = LE_force(1,j) * FiY(1,j);
        FLeZ (1,j) = LE_force(1,j) * FiZ(1,j);
    
   
    end


    % Moments

    for j=1:N-1
    
        MLeX (1,j)=  FLeY(j)*Z_C(1,j) - FLeZ(j)*Y_C(1,j);
        MLeY (1,j)= -FLeX(j)*Z_C(1,j) - FLeZ(j)*X_C(1,j);
        MLeZ (1,j)=  FLeX(j)*Y_C(1,j) - FLeY(j)*X_C(1,j);
    
    end

    MLeX = sum(MLeX);
    MLeY = sum(MLeY);
    MLeZ = sum(MLeZ);

    % Force

    FLeX = sum (FLeX);
    FLeY = sum (FLeY);
    FLeZ = sum (FLeZ);
    
end

Sail.FLeX=FLeX;
Sail.FLeY=FLeY;
Sail.FLeZ=FLeZ;
Sail.MLeX=MLeX;
Sail.MLeY=MLeY;
Sail.MLeZ=MLeZ;

