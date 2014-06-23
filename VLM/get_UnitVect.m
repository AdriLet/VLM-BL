function [Ex,Ey,Ez] = get_UnitVect (X1,Y1,Z1,X2,Y2,Z2)

Vx=X2-X1;
Vy=Y2-Y1;
Vz=Z2-Z1;

Ex= Vx/(Vx^2+Vy^2+Vz^2)^0.5;
Ey= Vy/(Vx^2+Vy^2+Vz^2)^0.5;
Ez= Vz/(Vx^2+Vy^2+Vz^2)^0.5;