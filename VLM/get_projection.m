function [Proj_UontoV,Proj_UontoV_x,Proj_UontoV_y,Proj_UontoV_z] = get_projection (Ux,Uy,Uz,Vx,Vy,Vz)

% Projects of U onto V. The output is a vector which is the component of U
% in the direction of V
% Rem: U and V are matrices.

N = max([length(Ux(1,:));length(Vx(1,:))]);
M = max([length(Ux(:,1));length(Vx(:,1))]);
    

norm_dot=((Ux.*Vx+Uy.*Vy+Uz.*Vz)./(Vx.*Vx+Vy.*Vy+Vz.*Vz));
 
Proj_UontoV_x= norm_dot.*Vx;
Proj_UontoV_y= norm_dot.*Vy;
Proj_UontoV_z= norm_dot.*Vz;

Proj_UontoV = ((Proj_UontoV_x).^2+(Proj_UontoV_y).^2+(Proj_UontoV_z).^2).^0.5;

for i=1:M
            
    for j=1:N
           
        if norm_dot(i,j) > 0    
                       
            Proj_UontoV (i,j)= Proj_UontoV(i,j);
                  
        end

        if  norm_dot (i,j)< 0
            
            Proj_UontoV (i,j)= -Proj_UontoV(i,j);
               
        end

    end

end