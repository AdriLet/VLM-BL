function [u,v,w] = get_V_from_semiinfvortl (x,y,z,x1,y1,z1,xE,yE,zE,gamma)

%Rem: E is a point on the semi-inf vortex line.

% clear all
% clc
% clf
% hold off
% close all
% 
% gamma=1;
% 
% x=0.83
% y=1
% z=0
% 
% x1=1
% y1=0.66
% z1=0
% 
% xE=10.9255
% yE=0.66
% zE=1.2187




%From my derivation

r1= ((x-x1)^2 + (y-y1)^2 + (z-z1)^2)^0.5;
einfx=(x1-xE)/(((x1-xE)^2+(y1-yE)^2+(z1-zE)^2)^0.5);
einfy=(y1-yE)/(((x1-xE)^2+(y1-yE)^2+(z1-zE)^2)^0.5);
einfz=(z1-zE)/(((x1-xE)^2+(y1-yE)^2+(z1-zE)^2)^0.5);

phix= (((y-y1)*einfz)-((z-z1)*einfy))   /  ((((y-y1)*einfz)-((z-z1)*einfy))^2+(((x-x1)*einfz)-((z-z1)*einfx))^2+(((x-x1)*einfy)-((y-y1)*einfx))^2);
phiy= -(((x-x1)*einfz)-((z-z1)*einfx))   /  ((((y-y1)*einfz)-((z-z1)*einfy))^2+(((x-x1)*einfz)-((z-z1)*einfx))^2+(((x-x1)*einfy)-((y-y1)*einfx))^2);
phiz= (((x-x1)*einfy)-((y-y1)*einfx))  /  ((((y-y1)*einfz)-((z-z1)*einfy))^2+(((x-x1)*einfz)-((z-z1)*einfx))^2+(((x-x1)*einfy)-((y-y1)*einfx))^2);

omega= -(((x-x1)*einfx+(y-y1)*einfy+(z-z1)*einfz)/(r1))+1;


u = gamma * phix * omega/(4*pi);
v = gamma * phiy * omega/(4*pi);
w = gamma * phiz * omega/(4*pi);

