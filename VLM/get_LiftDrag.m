function [L,Di] = get_LiftDrag (AWRef,FtotX,FtotY,FtotZ)

% Projection of the F_LE to have the lift and drag

% 1 The normal to the apparent wind

m=AWRef.Z /AWRef.X;
n=-1/m;
    
if n<0&&m~=0
    nAW_X_athref=-1/(1*1+n*n)^0.5;
    nAW_Z_athref=-n/(1*1+n*n)^0.5;
end

if n>0&&m~=0
    nAW_X_athref=1/(1*1+n*n)^0.5;
    nAW_Z_athref=n/(1*1+n*n)^0.5;   
end

if m==0
    nAW_X_athref=0;                          
    nAW_Z_athref=1;          
end
      

% 2 The projection

[L] = get_projection (FtotX,FtotY,FtotZ,   nAW_X_athref, 0, nAW_Z_athref);
[Di] = get_projection (FtotX,FtotY,FtotZ,  AWRef.X,0,AWRef.Z);