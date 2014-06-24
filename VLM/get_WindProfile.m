function [WindProfile, AWRef]=get_WindProfile(TWSRef,TWA,Vs,Lambda,Windshear)
% This function compute the wind profile at the same height than Flow

% 1. Height from FLow
h=[1;3.33333300000000;5;10;15;20;25;30;35;40;45;50;55;60;65;70;75;];

% 2 The roughness length
y0 = 0.00005 * TWSRef^2 /9.81;

% 3 The Apparent Wind (From Hansen's Tesis)
if Windshear==1
    TWS(:,1) =  TWSRef * (log(h./y0))/(log(10/y0));
end

if Windshear==0
    TWS(:,1)=TWSRef*ones(length(h),1);
end

AW_X(:,1)=  TWS(:,1) * cos (TWA*pi/180) + Vs*cos(-Lambda*pi/180);
AW_Z(:,1)=  TWS(:,1) * sin (TWA*pi/180)+ Vs*sin(-Lambda*pi/180);

AWA=atan(AW_Z./AW_X)*180/pi;
AWS=(AW_X.^2+AW_Z.^2).^0.5;

% 4 The Apparent Wind at 10m (=heigth reference).
AWRef_X= TWSRef(:,1) * cos (TWA*pi/180) + Vs*cos(-Lambda*pi/180);
AWRef_Z= TWSRef(:,1) * sin (TWA*pi/180) + Vs*sin(-Lambda*pi/180);

WindProfile=[h,TWS,TWA*ones(length(h),1),AWA,AWS];
AWRef.X=AWRef_X;
AWRef.Z=AWRef_Z;
AWRef.Mgn=(AWRef_X^2 + AWRef_Z^2)^0.5;

end
