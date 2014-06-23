function CurrentLine= getStreamLine2 (Sail, h, in)
%options
interMethod='nearest';



% This function follow the stream line at height = h. in is 1 if we are
% looking at intrado, 0 if looking at extrado.

X_C=Sail.X_C;
Y_C=Sail.Y_C;
Z_C=Sail.Z_C;
X_V=Sail.X_V;
Y_V=Sail.Y_V;
Z_V=Sail.Z_V; 
GAMMA=Sail.GAMMA;

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


  i=min(find( Y_C(1,:) > h ) ); %élement juste supérieur à h dans la première ligne horizontale (guindant).
  if isempty(i)
        CurrentLine.X=0;
        CurrentLine.Y=0;
        CurrentLine.Z=0;
        CurrentLine.V=0;
        return

  elseif  i==1 
    CurrentLine.X=0;
    CurrentLine.Y=0;
    CurrentLine.Z=0;
    CurrentLine.V=0;
    return
  end
  %j=min(find( Vx(:,i)>0  ))+2;
  pente_etais=(Y_C(1,i)-Y_C(1,i-1))/(X_C(1,i)-X_C(1,i-1));
    lambda=(h-Y_C(1,i))/(Y_C(1,i-1)-Y_C(1,i));
  j=min(find(      ((lambda*Vy(:,i-1)+(1-lambda)*Vy(:,i))   ./   (lambda*Vx(:,i-1)+(1-lambda)*Vx(:,i)) < pente_etais ) & (Vx(:,i)>0)  ));

    lambda=(h-Y_C(j,i))/(Y_C(j,i-1)-Y_C(j,i));

  Xfl=zeros(1,4*M-1);
  Yfl=zeros(1,4*M-1);
  Zfl=zeros(1,4*M-1);
  Vfl=zeros(1,4*M-1);
Xfl(1,1)=lambda*X_C(j,i-1)+(1-lambda)*X_C(j,i);
Yfl(1,1)=lambda*Y_C(j,i-1)+(1-lambda)*Y_C(j,i);
Zfl(1,1)=lambda*Z_C(j,i-1)+(1-lambda)*Z_C(j,i);
Vfl(1,1)=lambda*V(j,i-1)+(1-lambda)*V(j,i);
%V0=[lambda*Vx(1,i-1)+(1-lambda)*Vx(1,i);  lambda*Vy(1,i-1)+(1-lambda)*Vy(1,i); lambda*Vz(1,i-1)+(1-lambda)*Vz(1,i)];
  V0= [griddata(X_C,Y_C,Vx,Xfl(1,1),Yfl(1,1),interMethod); griddata(X_C,Y_C,Vy,Xfl(1,1),Yfl(1,1),interMethod) ; griddata(X_C,Y_C,Vz,Xfl(1,1),Yfl(1,1),interMethod)];

    
    %3.2 delta x
    
    deltas=10*abs(X_C(1,i)-X_C(2,i));
    
    %%%%%%
   j=1;

while and(isOnTheShape(Xfl(1,j),Yfl(1,j),Zfl(1,j),Sail,deltas),j<4*M)
    j=j+1;

    x= [ Xfl(1,j-1);  Yfl(1,j-1);  Zfl(1,j-1)]+deltas*V0/sqrt(sum(V0.^2));
        Xfl(1, j)=x(1);
        Yfl(1, j)=x(2);
        Zfl(1, j) =x(3);
    Zfl(1,j)=griddata(X_C,Y_C,Z_C,Xfl(1,j),Yfl(1,j),'linear');
     V0= [griddata(X_C,Y_C,Vx,Xfl(1,j),Yfl(1,j),interMethod); griddata(X_C,Y_C,Vy,Xfl(1,j),Yfl(1,j),interMethod) ; griddata(X_C,Y_C,Vz,Xfl(1,j),Yfl(1,j),interMethod)];
    Vfl(1,j)=sum(V0.^2);

end







i=min(find(Yfl==0));
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







end











%% Local functions



    function bool=isOnTheShape(X,Y,Z,Sail,deltas)
    

  i=min(find( Sail.Y_C(end,:) > Y) );   
  if isempty(i) | i==1 | i==size(Sail.Y_C,2)
      bool=false;
      return
  end
          
    lambda=(Y-Sail.Y_C(end,i))/(Sail.Y_C(end,i-1)-Sail.Y_C(end,i))

               bool= X<lambda*Sail.X_C(end,i-1)+(1-lambda)*Sail.X_C(end,i)-deltas;                       %if this is true, the points get ou by luff
          

       bool=bool &max(max( Sail.X_C))>X & min(min(Sail.X_C))<X & max(max(Sail.Y_C))>Y & min(min(Sail.Y_C))<Y;
 
        
    end
%end