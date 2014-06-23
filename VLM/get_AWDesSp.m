function [AWA,AWS]=get_AWDesSp(Y,TWSRef,TWA,Vs,Lambda,Windshear)
% 1 The roughness length
y0 = 0.00005 * TWSRef^2 /9.81;

% 2 The Apparent Wind (From Hansen's Tesis)
if Windshear==1
    TWS(:,:) =  TWSRef * (log(Y(1,:)./y0))/(log(10/y0));
end

if Windshear==0
    TWS(:,:)=TWSRef*ones(M-1,N-1);
end

% Note that AW_X is the apparent wind at every horizontal panel
AW_X(1,:)=  TWS .* cos (TWA*pi/180) + Vs*cos(-Lambda*pi/180);
AW_Z(1,:)=  TWS .* sin (TWA*pi/180)+ Vs*sin(-Lambda*pi/180);

AWA=atan(AW_Z./AW_X)*180/pi;
AWS=(AW_X.^2+AW_Z.^2).^0.5;