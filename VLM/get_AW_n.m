function [Sail] = get_AW_n (Sail,TWSRef,TWA,Vs,Lambda,Windshear)
% This function projects the apparent wind on the normal of the sail panels
M =length(Sail.X_V(:,1));
N =length(Sail.X_V(1,:));

% 2 The roughness length
y0 = 0.00005 * TWSRef^2 /9.81;

% 3 The Apparent Wind (From Hansen's Tesis)
if Windshear==1
    TWS(:,:) =  TWSRef * (log(Sail.Y_C./y0))/(log(10/y0));
end

if Windshear==0
    TWS(:,:)=TWSRef*ones(M-1,N-1);
end

% Note that AW_X is the apparent wind at every panel of the Sail
AW_X(:,:)=  TWS .* cos (TWA*pi/180) + Vs*cos(-Lambda*pi/180);
AW_Z(:,:)=  TWS .* sin (TWA*pi/180)+ Vs*sin(-Lambda*pi/180);
AW_Y(:,:)=zeros(M-1,N-1);

[AW_n] = get_projection (AW_X,AW_Y,AW_Z,Sail.nx,Sail.ny,Sail.nz);

Sail.AW_n=AW_n;
Sail.AW_X=AW_X;
Sail.AW_Z=AW_Z;
