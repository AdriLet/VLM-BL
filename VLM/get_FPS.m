function [Sail] = get_FPS(Sail,AWRef,rho)

X_V=Sail.X_V;
Y_V=Sail.Y_V;
Z_V=Sail.Z_V;

X_C=Sail.X_C;
Y_C=Sail.Y_C;
Z_C=Sail.Z_C;
GAMMA=Sail.GAMMA;

Vmeanx=Sail.Vmeanx;
Vmeany=Sail.Vmeany;
Vmeanz=Sail.Vmeanz;

nx=Sail.nx;
ny=Sail.ny;
nz=Sail.nz;

M=length(X_V(:,1));
N=length(X_V(1,:));

% 1 Gradient of the Vortex Strength (doublicity)
% 1.1 s and t coordinates
% 1.1.1 Origin at the tack of the Sail
X_V0=X_V-X_V(1,1);
Y_V0=Y_V-Y_V(1,1);
Z_V0=Z_V-Z_V(1,1);

X_C0=X_C-X_V(1,1);
Y_C0=Y_C-Y_V(1,1);
Z_C0=Z_C-Z_V(1,1);

% 1.1.2 Unheeled the datas (rotation around x)
theta = atan (Z_V0(1,end)/Y_V0(1,end));

Y_V0_rot =  Y_V0 * cos(theta) + Z_V0 * sin (theta);
Z_V0_rot =  - Y_V0 * sin(theta) + Z_V0 * cos (theta);

Y_C0_rot =  Y_C0 * cos(theta) + Z_C0 * sin (theta);
Z_C0_rot =  - Y_C0 * sin(theta) + Z_C0 * cos (theta);

%1.1.3 Origin at the tack  + Adim (X and Y_rot) 
for j=1:N-1    
    s_C (:,j) = (X_C0(:,j)-X_C0(1,j)) / (X_C0(end,j)-X_C0(1,j));                             
end

for i=1:M-1
    t_C (i,:) = (Y_C0_rot(i,:)-Y_C0_rot(i,1)) / (Y_C0_rot(i,end)-Y_C0_rot(i,1));                          
end

% 1.2 The gradient in 3D
%1.2.1 Computations
    %1 dGAMMA/ds
    %-----------
    [dGAMMAds]= get_grad_alongchord(s_C,t_C,GAMMA);

    %2 dGAMMA/dt
    %-----------
    [dGAMMAdt]= get_grad_alongspan(s_C,t_C,GAMMA);

    %3 d/ds
    %------
    %3.1 dx/ds
    [dxds]= get_grad_alongchord(s_C,t_C,X_C);

    %3.2 dy/ds
    [dyds]= get_grad_alongchord(s_C,t_C,Y_C);

    %3.3 dz/ds
    [dzds]= get_grad_alongchord(s_C,t_C,Z_C);

    %4 d/dt
    %------
    %4.1 dx/dt
    [dxdt]= get_grad_alongspan(s_C,t_C,X_C);

    %4.2 dy/dt
    [dydt]= get_grad_alongspan(s_C,t_C,Y_C);

    %4.3 dz/dt
    [dzdt]= get_grad_alongspan(s_C,t_C,Z_C);

    %4 d/du
    %------
    dxdu=nx;
    dydu=ny;
    dzdu=nz;
    
%1.2.2 Solve the system
for i=1:M-1    
    for j=1:N-1
        A=[dxds(i,j) dyds(i,j) dzds(i,j);...
           dxdt(i,j) dydt(i,j) dzdt(i,j);...
           dxdu(i,j) dydu(i,j) dzdu(i,j)];

        B=[dGAMMAds(i,j);dGAMMAdt(i,j);0];

        C=A\B;
        
        dGAMMAdx(i,j)=C(1);
        dGAMMAdy(i,j)=C(2);
        dGAMMAdz(i,j)=C(3);        
    end    
end


% 2 The Velocities and Cps
% 2.1 Velocity upper
Vxup = Vmeanx - 0.5 * dGAMMAdx;
Vyup = Vmeany - 0.5 * dGAMMAdy;
Vzup = Vmeanz - 0.5 * dGAMMAdz;

for i=1:M-1    
    for j=1:N-1        
        Vup (i,j) = (Vxup(i,j)^2 + Vyup(i,j)^2 + Vzup(i,j)^2)^0.5;
        Cpup (i,j) = 1 - (Vup (i,j)/AWRef.Mgn)^2;       
    end    
end

% 2.2 Velocity lower
Vxlo = Vmeanx + 0.5 * dGAMMAdx;
Vylo = Vmeany + 0.5 * dGAMMAdy;
Vzlo = Vmeanz + 0.5 * dGAMMAdz;

for i=1:M-1   
    for j=1:N-1       
        Vlo (i,j) = (Vxlo(i,j)^2 + Vylo(i,j)^2 + Vzlo(i,j)^2)^0.5;
        Cplo (i,j) = 1 - (Vlo (i,j)/AWRef.Mgn)^2;        
    end    
end

% 3 The Pressure Jump
for i=1:M-1    
    for j=1:N-1        
        dp(i,j) = 0.5 * rho * ( Vup(i,j)^2 - Vlo(i,j)^2);        
    end    
end
   
% 4 The Cp
dCp= Cplo-Cpup;

% 5 The force and the moment
for j=1:N-1    
    for i=1:M-1      
        FPsX (i,j) = Sail.area(i,j) * dp(i,j)* nx(i,j);
        FPsY (i,j) = Sail.area(i,j) * dp(i,j)* ny(i,j);
        FPsZ (i,j) = Sail.area(i,j) * dp(i,j)* nz(i,j);
        
        MPsX (i,j) =   Y_C(i,j) * FPsZ (i,j) - Z_C(i,j) * FPsY (i,j);
        MPsY (i,j) = - X_C(i,j) * FPsZ (i,j) + Z_C(i,j) * FPsX (i,j);
        MPsZ (i,j) =   X_C(i,j) * FPsY (i,j) - Y_C(i,j) * FPsX (i,j);                            
    end    
end

FPsX=sum(sum(FPsX));
FPsY=sum(sum(FPsY));
FPsZ=sum(sum(FPsZ));

MPsX=sum(sum(MPsX));
MPsY=sum(sum(MPsY));
MPsZ=sum(sum(MPsZ));

Sail.FX=FPsX;
Sail.FY=FPsY;
Sail.FZ=FPsZ;
Sail.MX=MPsX;
Sail.MY=MPsY;
Sail.MZ=MPsZ;
Sail.dCp_C=dCp;
Sail.Cplo_C=Cplo;
Sail.Cpup_C=Cpup;

%     figure (4)
%     surface(X_C,Y_C,dCp); view(-36,36);drawnow
%     title('dCp');
% %     axis('equal')
%     xlabel('X [m]')
%     ylabel('Y [m]')
%     zlabel('dCp [-] ')
