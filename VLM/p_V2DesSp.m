function [dCp_DesSp]=p_V2DesSp(sail,XDesSp,YDesSp,ZDesSp)
% the desired spacing has to be another vortex lines discretisation
% 1 Generation of the s and t coordinate
[sail.s_V,sail.t_V] = get_sandt (sail.X_V,sail.Y_V,sail.Z_V);
[sDesSp,tDesSp] = get_sandt (XDesSp,YDesSp,ZDesSp);

% 2 Interpolation of the pressure in s and t coordinates
dCp_DesSp = griddata(sail.s_V,sail.t_V,sail.dCp_V,sDesSp,tDesSp,'cubic');
