function CurrentLine= getStreamLine4(Sail, j, in)
%options
interMethod='natural';
j1=j;


% This function follow the stream line at height = h. in is 1 if we are
% looking at intrado, 0 if looking at extrado.

X_C=Sail.X_C;
Y_C=Sail.Y_C;
Z_C=Sail.Z_C;
X_V=Sail.X_V;
Y_V=Sail.Y_V;
Z_V=Sail.Z_V; 
GAMMA=Sail.GAMMA;
X=Sail.X;
Y=Sail.Y;
Z=Sail.Z;

Vmeanx=Sail.Vmeanx;
Vmeany=Sail.Vmeany;
Vmeanz=Sail.Vmeanz;

nx=Sail.nx;
ny=Sail.ny;
nz=Sail.nz;

M=length(X_V(:,1));
N=length(X_V(1,:));


%% 1 Gradient of the Vortex Strength (doublicity)
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


%% 2 The Velocities 
% 2.1 Velocity on In/Out of the sail
Vx = Vmeanx +(-1+2*in)*0.5 * dGAMMAdx;        % if in is 1 then intrado then Vmean+0.5GAMMA
Vy = Vmeany +(-1+2*in)*0.5 * dGAMMAdy;        % elsein=0 and Vmean-0.5GAMMA
Vz = Vmeanz +(-1+2*in)*0.5 * dGAMMAdz;  

V  = (Vx.^2 + Vy.^2 + Vz.^2).^0.5;






                

%% 3 The line
% 3.1 Starting point 


j=j1;
  %j=min(find( Vx(:,i)>0  ))+2;
  pente_etais=(Y(1,1)-Y(1,2)) /  (X(1,1)-X(1,2))-70*pi/180;

 i=min( find(  (Vy(:,j)  ./   Vx(:,j) < pente_etais ) & (Vx(:,j)>0)  ));


  Xfl=zeros(1,4*M-1);
  Yfl=zeros(1,4*M-1);
  Zfl=zeros(1,4*M-1);
  Vfl=zeros(1,4*M-1);
Xfl(1,1)=X_C(i,j);
Yfl(1,1)=Y_C(i,j);
Zfl(1,1)=Z_C(i,j);
Vfl(1,1)=V(i,j);
%V0=[lambda*Vx(1,i-1)+(1-lambda)*Vx(1,i);  lambda*Vy(1,i-1)+(1-lambda)*Vy(1,i); lambda*Vz(1,i-1)+(1-lambda)*Vz(1,i)];
  V0= [Vx(i,j) ; Vy(i,j) ; Vz(i,j)];

    
    %3.2 delta x
    
    deltas=0.2*abs(X_V(3,j)-X_V(2,j));
    
    %%%%%%
   j=1;

while and(isOnTheShape(Xfl(1,j),Yfl(1,j),Zfl(1,j),Sail,deltas),j<15*M)
    j=j+1;

    x= [ Xfl(1,j-1);  Yfl(1,j-1);  Zfl(1,j-1)]+deltas*V0/sqrt(sum(V0.^2));
        Xfl(1, j)=x(1);
        Yfl(1, j)=x(2);
        Zfl(1, j) =x(3);
    Zfl(1,j)=griddata(X_C,Y_C,Z_C,Xfl(1,j),Yfl(1,j),'linear');
     V0= [griddata(X_C,Y_C,Vx,Xfl(1,j),Yfl(1,j),interMethod); griddata(X_C,Y_C,Vy,Xfl(1,j),Yfl(1,j),interMethod) ; griddata(X_C,Y_C,Vz,Xfl(1,j),Yfl(1,j),interMethod)];
    Vfl(1,j)=sqrt(sum(V0.^2));

end










i=min(find(Yfl==0 | isnan(Zfl) ));
if isempty(i)
    CurrentLine.X=Xfl;
    CurrentLine.Y=Yfl;
    CurrentLine.Z=Zfl;
    CurrentLine.V=Vfl;
else
    CurrentLine.X=Xfl(1:i-1);
    CurrentLine.Y=Yfl(1:i-1);
    CurrentLine.Z=Zfl(1:i-1);
    CurrentLine.V=Vfl(1:i-1);
end

h=mean(Yfl(1,:));
CurrentLine.Vinf=spline(Sail.outputs.WindProfile.Height,Sail.outputs.WindProfile.AWS,h);




end











%% Local functions



    function bool=isOnTheShape(X,Y,Z,Sail,deltas)
    

 j=min(find( Sail.Y(end,:) > Y) );   
  if isempty(j) | j==1 | j==size(Sail.Y,1)
      bool=false;
      return
  end
          
    lambda=(Y-Sail.Y(end,j))/(Sail.Y(end,j-1)-Sail.Y(end,j));

               bool= X<lambda*Sail.X(end,j-1)+(1-lambda)*Sail.X(end,j)-deltas;                       %if this is true, the points get ou by luff
          

       bool=bool &max(max( Sail.X))>X & min(min(Sail.X))<X & max(max(Sail.Y))>Y & min(min(Sail.Y))<Y;
 
        
    end
%end