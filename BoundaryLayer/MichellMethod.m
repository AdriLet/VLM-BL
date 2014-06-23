function [stall,line]=MichellMethod(Streamline)
% This function takes a streamline with X, Y, Z, V in it. Then it calculates theta, Cf, delta, delta1 and put it in line. Stall is true when stall is detected on the stream line. It occurs when the form factor  H=3.5 The turbulent boundary layer method is a 2D integral method named Michell et al 's method (1969)


%constantes
nu=15.6*10^(-6);           %dynamic viscosity
chi=0.41;                     %Von Karman's Constant
relaxFactor=0.4;


Vinf=Streamline.Vinf;
N=max(size(Streamline.X));
theta=zeros(1,N);
delta1=zeros(1,N);
delta=zeros(1,N);
H=zeros(1,N);
Cf=zeros(1,N);
alpha=zeros(1,N);

X=Streamline.X;
Y=Streamline.Y;
Z=Streamline.Z;
V=Streamline.V;
dS(2:N)=sqrt((X(1:N-1)-X(2:N)).^2+(Y(1:N-1)-Y(2:N)).^2+(Z(1:N-1)-Z(2:N)).^2);          % S is curvilinear coordinate.
dS(1)=0;

alphaTest=2;
RS=Vinf*dS(2)/nu;

delta(2)=0.475*dS(2)/RS^(1/5);
theta(2)=7/72*delta(2);
delta1(2)=delta(2)/8;

deltaMinusDelta1=delta(2)-delta1(2);
Theta=theta(2);



 syms x positive
for i=2:N

    while    abs(alpha(i)-alphaTest)>10^(-3)

        %syms x positive
        H(i) =min(eval( solve(  deltaMinusDelta1/Theta== x*(alphaTest*x+1)/(x-1)  , x)));
        delta1(i)=H(i)*Theta;
       delta(i)=deltaMinusDelta1+delta1(i);
       %sysm x positive
       Rdelta1=V(i)*delta1(i)/nu;
       G=eval( solve( x*H(i)/(H(i)-1) ==1/chi*log(Rdelta1)+1.4*x+20/x-8.33 , x));
       F1=((0.613*G-(3.6+76.86*(1/G-0.154)^2)/G))/(1+10/G^2);
        alpha(i)=alphaTest;
       alphaTest=(G/F1-1)*relaxFactor+alpha(i)*(1-relaxFactor);
        gamma=(H(i)-1)/(H(i)*G);
        Cf(i)=2*gamma^2;
        
        
    end
DVeDs=(V(i)-V(i-1))/dS(i);
    
DthetaDs=Cf(i)/2-theta(i)*(H(i)+2)/V(i)*DVeDs;
P=0.558*G^2/F1-1.325/F1;
CE=gamma*(0.074*G-1.0957/G);
DdeltaMinusDelta1Ds=CE-deltaMinusDelta1/V(i)*DVeDs;
   deltaMinusDelta1 =deltaMinusDelta1 +DdeltaMinusDelta1Ds *dS(i);
   theta(i)=Theta;
    Theta=Theta+DthetaDs*dS(i);
    
end







end
