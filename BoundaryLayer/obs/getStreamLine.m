function CurrentLine= getStreamLine (Sail, h, in)
% This function follow the flow line at height = h. in is 1 if we are
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

dx=X_C(1,end)-X_C(2,end);

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
  if i==1 
    CurrentLine.X=0;
    CurrentLine.Y=0;
    CurrentLine.Z=0;
    CurrentLine.V=0;
    return
  end
  lambda=(h-Y_C(1,i))/(Y_C(1,i-1)-Y_C(1,i));
  Xfl=zeros(1,M-1);
  Yfl=zeros(1,M-1);
  Zfl=zeros(1,M-1);
  Vfl=zeros(1,M-1);
Xfl(1,1)=lambda*X_C(1,i-1)+(1-lambda)*X_C(1,i);
Yfl(1,1)=lambda*Y_C(1,i-1)+(1-lambda)*Y_C(1,i);
Zfl(1,1)=lambda*Z_C(1,i-1)+(1-lambda)*Z_C(1,i);
Vfl(1,1)=lambda*V(1,i-1)+(1-lambda)*V(1,i);
%V0=[lambda*Vx(1,i-1)+(1-lambda)*Vx(1,i);  lambda*Vy(1,i-1)+(1-lambda)*Vy(1,i); lambda*Vz(1,i-1)+(1-lambda)*Vz(1,i)];
    V0= [griddata(X_C,Y_C,Vx,Xfl(1,1),Yfl(1,1),'nearest'); griddata(X_C,Y_C,Vy,Xfl(1,1),Yfl(1,1),'nearest') ; griddata(X_C,Y_C,Vz,Xfl(1,1),Yfl(1,1),'nearest')];

    
    %3.2 delta x
    
    deltax=10*max(X_C)-min(X_C);
   

for j=2:M-1
    
    [Xfl(1,j), Yfl(1,j), Zfl(1,j)]=getClosestFromTheLine([X_C(j,i-1); Y_C(j,i-1); Z_C(j,i-1)], [X_C(j,i); Y_C(j,i); Z_C(j,i)], V0, [Xfl(1); Yfl(j-1); Zfl(j-1)]);
     i= min(find( Y_C(j,:) > Yfl(1,j) ) );
     if i==1
         %disp 'break1'
         break
     elseif isempty(i)==1 
        %disp 'break0'
         break
     end
     
     lambda= (Yfl(1,j)-Y_C(j,i))/(Y_C(j,i-1)-Y_C(j,i));
     Xfl(j)=[lambda*X_C(j,i-1)+(1-lambda)*X_C(j,i)];
     Yfl(j)=[lambda*Y_C(j,i-1)+(1-lambda)*Y_C(j,i)];
     Zfl(j)=[lambda*Z_C(j,i-1)+(1-lambda)*Z_C(j,i)];
   Vfl(1,j)= lambda*V(j,i-1)+(1-lambda)*V(j,i) ;
    %V0= [lambda*Vx(j,i-1)+(1-lambda)*Vx(j,i);  lambda*Vy(j,i-1)+(1-lambda)*Vy(j,i); lambda*Vz(j,i-1)+(1-lambda)*Vz(j,i)];
     %V0=[interp3(X_C,Y_C,Z_C,Vx,Xfl(1,j),Yfl(1,j),Zfl(1,j)); interp3(X_C,Y_C,Z_C,Vy,Xfl(1,j),Yfl(1,j),Zfl(1,j)); interp3(X_C,Y_C,Z_C,Vz,Xfl(1,j),Yfl(1,j),Zfl(1,j))];
    V0= [griddata(X_C,Y_C,Vx,Xfl(1,j),Yfl(1,j),'nearest'); griddata(X_C,Y_C,Vy,Xfl(1,j),Yfl(1,j),'nearest') ; griddata(X_C,Y_C,Vz,Xfl(1,j),Yfl(1,j),'nearest')];
%  
%    if or( X_C(end,j) < Xfl(1,j), Xfl(1,j) < X_C(1,j) )
%        Xfl(1,j)=0;
%        Yfl(1,j)=0;
%        Zfl(1,j)=0;
%        Vfl(1,j)=0;
%    elseif or(Y_C(1,j)>Yfl(1,j),Yfl(1,j)>Y_C(end,j))
%       Xfl(1,j)=0;
%        Yfl(1,j)=0;
%        Zfl(1,j)=0;
%        Vfl(1,j)=0;
%      elseif or(Z_C(1,j)>Zfl(1,j),Zfl(1,j)>Z_C(end,j))  
%        Xfl(1,j)=0;
%        Yfl(1,j)=0;
%        Zfl(1,j)=0;
%        Vfl(1,j)=0;     
%    end
    
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
% This function finds the closest point from line (Xi1Xi2), on the line
% directed by V0, starting in X0. It returns Xout, the coordinates of this
% point. X is always a (3,1) vector
function [Xout, Yout, Zout]=getClosestFromTheLine (Xi1, Xi2, V0, X0)



T=[cross(Xi2-Xi1,V0) (Xi2-Xi1)  V0]\(X0-Xi2);
t2=T(2);
out=t2*(Xi2-Xi1)+Xi2;
Xout=out(1);
Yout=out(2);
Zout=out(3);
end

%end