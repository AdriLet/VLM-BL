function [u,v,w] = get_V_from_vortl (x,y,z,x1,y1,z1,x2,y2,z2,gamma)

% clear all
% clc
% clf
% hold off
% close all
% 
% 
% x=0.16
% y=0.33
% z=0.055
% 
% x1=0
% y1=0
% z1=0
% 
% x2=0.33
% y2=0
% z2=0.0888
% 
% gamma=1


  
 %From Mason
 a=(y-y1)*(z-z2)-(z-z1)*(y-y2);
 b=(x-x1)*(z-z2)-(z-z1)*(x-x2);
 c=(x-x1)*(y-y2)-(y-y1)*(x-x2);
 denum=a^2+b^2+c^2;

 phix=a/denum;
 phiy=-b/denum;
 phiz=c/denum;
 
 omega=(((x2-x1)*(x-x1)+(y2-y1)*(y-y1)+(z2-z1)*(z-z1))/(((x-x1)^2+(y-y1)^2+(z-z1)^2)^0.5))-(((x2-x1)*(x-x2)+(y2-y1)*(y-y2)+(z2-z1)*(z-z2))/(((x-x2)^2+(y-y2)^2+(z-z2)^2)^0.5));
 
 u=gamma*phix*omega / (4*pi);
 v=gamma*phiy*omega / (4*pi);
 w=gamma*phiz*omega / (4*pi);