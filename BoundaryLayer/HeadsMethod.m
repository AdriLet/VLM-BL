function [ stall, Xstall, Ystall, Zstall, Cf, H, theta] = HeadsMethod( StreamLine )
%It calculates boundary layer according to Head's method, based on Ludwieg-Tillman Cf
%if stall is detected, stall is true and X,Y,Z contains the point on wich
%stall apears.
nu=15.6*10^(-6);           %dynamic viscosity

X=StreamLine.X;
Y=StreamLine.Y;
Z=StreamLine.Z;
V=StreamLine.V;
Vinf=StreamLine.Vinf;
N=max(size(V));
if N==1
    
    stall=false;
    Xstall=0;
    Ystall=0;
    Zstall=0;
    Cf=zeros(1,N);
    H=zeros(1,N);
    theta=zeros(1,N);
    return
end

dS(2:N)=sqrt((X(1:N-1)-X(2:N)).^2+(Y(1:N-1)-Y(2:N)).^2+(Z(1:N-1)-Z(2:N)).^2);          % S is curvilinear coordinate.
dS(1)=0;
S(1)=0;
for i=2:N
    S(i)=S(i-1)+dS(i);
end


theta=zeros(1,N);
delta1=zeros(1,N);
delta=zeros(1,N);


H=zeros(1,N);
Cf=zeros(1,N);
RS=Vinf*dS(2)/nu;
delta(1)=0.475*dS(2)/RS^(1/5);
theta(1)=7/72*delta(1);
delta1(1)=delta(1)/8;

H(1)=delta1(1)/theta(1);
dVds(2:N)=(V(2:N)-V(1:N-1))./dS(2:N);
dVdS(1)=0;

dThetadS=0;
dHdS=0;
R_theta=Vinf*theta(1)/nu;

    Cf(1)=0.246*10^(-0.678*H(1))*R_theta^(-0.268);


for j=2:N

    H_star=1.535*(H(j-1)-0.7)^(-2.715)+3.3;
    CE=0.0306*(H_star-3)^(-0.653);    
    dThetadS=Cf(j-1)/2-theta(j-1)*(H(j-1)+2)/V(j-1)*dVds(j);
    dHdS=1/(theta(j-1)*1.535*(H(j-1)-0.7)^(-3.715)*(-2.715))*(CE-theta(j-1)*H_star/V(j-1)*dVds(j)-H_star*dThetadS);
        H(j)=H(j-1)+dS(j)*dHdS;
    theta(j)=theta(j-1)+dS(j)*dThetadS;
     R_theta=Vinf*theta(j)/nu;
    Cf(j)=0.246*10^(-0.678*H(j))*R_theta^(-0.268);
end

i=min(find(H>2.35));
if isempty(i)
    stall=false;
    Xstall=nan;
    Ystall=nan;
    Zstall=nan;
else
    stall=true;
    Xstall=X(i);
    Ystall=Y(i);
    Zstall=Z(i);
end



end



