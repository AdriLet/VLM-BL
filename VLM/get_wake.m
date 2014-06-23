function [Sail] = get_wake (Sail,TWSRef,TWA,Vs,Lambda,Windshear)
% 4-3-2011: Add the windshear option

M =length(Sail.X_V(:,1));
N =length(Sail.X_V(1,:));

L_wake=5;                                                                 %=20 By default, but is for the display because in the code is infinity
WAKE_X=zeros(10*M,N);
WAKE_Y=zeros(10*M,N);
WAKE_Z=zeros(10*M,N);

% 1 Creation of Wake
for j=1:N
    WAKE_X(:,j)=linspace(0,L_wake,10*M);
end

for i=1:10*M
    WAKE_Y(i,:)=Sail.Y_V(end,:);
end

% 2 Rotation in the apparent wind of the trealing edge cell
% 2.1 The roughness length
y0 = 0.00005 * TWSRef^2 /9.81;

% 2.2 The Apparent Wind (From Hansen's Tesis)
if Windshear==1
    TWS(1,:) =  TWSRef * (log(Sail.Y_V(end,:)./y0))/(log(10/y0));
end

if Windshear==0
    TWS(1,:)=TWSRef*ones(1,:);
end

% Note that AW_X is the apparent wind at every horizontal panel
AW_X(1,:)=  TWS .* cos (TWA*pi/180) + Vs*cos(-Lambda*pi/180);
AW_Z(1,:)=  TWS .* sin (TWA*pi/180)+ Vs*sin(-Lambda*pi/180);

% 3 Computation of the apparent wind angle at the trailing edge
theta=atan(AW_Z./AW_X);

% 4 Rotation
WAKE_X_rot=zeros(10*M,N);
WAKE_Z_rot=zeros(10*M,N);
for j=1:N
    WAKE_X_rot(:,j) = WAKE_X(:,j) * cos (theta(1,j)) - WAKE_Z (:,j)* sin (theta(1,j));
    WAKE_Z_rot(:,j) =  WAKE_X (:,j)* sin (theta(1,j)) + WAKE_Z(:,j) * cos (theta(1,j));
end

% 5 On the trailing edge
WAKE_X_rot_atTE=zeros(10*M,N);
WAKE_Z_rot_atTE=zeros(10*M,N);

for j=1:N
    WAKE_X_rot_atTE (:,j)= WAKE_X_rot (:,j) + Sail.X_V(end,j);
    WAKE_Z_rot_atTE (:,j)= WAKE_Z_rot (:,j) + Sail.Z_V(end,j);
end

Sail.wake_X=WAKE_X_rot_atTE;
Sail.wake_Y=WAKE_Y;
Sail.wake_Z=WAKE_Z_rot_atTE;