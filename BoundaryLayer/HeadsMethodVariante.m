function [ stall, Xstall, Ystall, Zstall] = HeadsMethodVariante( StreamLine )
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
Cf=zeros(1,N);
H1=zeros(1,N);


H=zeros(1,N);
RS=Vinf*dS(2)/nu;
delta(1)=0.475*dS(2)/RS^(1/5);
theta(1)=7/72*delta(1);
delta1(1)=delta(1)/8;

H(1)=delta1(1)/theta(1);
H1(1)=2+1.5*(1.12/(H(1)-1))^(1/0.915)+0.5*((H(1)-1)/1.12)^(1/915);
dVdS(2:N)=(V(2:N)-V(1:N-1))./dS(2:N);
dVdS(1)=0;

dThetadS=0;
dH1dS=0;



for j=2:N
    H1(j)=H1(j-1)+dS(j)*dH1dS;
    theta(j)=theta(j-1)+dS(j)*dThetadS;
    CE=0.0299*(H1(j)-3)^(-0.6169);
    R_theta=Vinf*theta(j)/nu;
    Cf0=0.012/(log(R_theta)/log(10)-0.64)-0.00093;
    H0=(1-6.8*sqrt(Cf0/2))^(-1);
    H(j)=1+1.12*(H1(j)-2-sqrt((H1(j)-2)^2-3))^0.915;
    Cf(j)=(0.9/(H(j)/H0-0.4)-0.5)*Cf0;
    
    dThetadS=Cf(j)/2-theta(j)*(H(j)+2)/V(j)*dVdS(j);
    dH1dS=1/theta(j)*(CE-H1(j)*(Cf(j)/2-(H(j)+1)*theta(j)/V(j)*dVdS(j)));
end

i=min(find(H>2.35));
if isempty(i)
    stall=false;
    Xstall=0;
    Ystall=0;
    Zstall=0;
else
    stall=true;
    Xstall=X(i);
    Ystall=Y(i);
    Zstall=Z(i);
end



end

