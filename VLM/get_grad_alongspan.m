function [dphidt]= get_grad_alongspan(s_C,t_C,phi)

M=length(s_C(:,1))+1;
N=length(s_C(1,:))+1;

% dphi/dt
%--------

for i=1:M-1
    
    for j=1:N-1      
        
        if j==1                                                            % grad @ the Foot
            
        dx1= t_C(i,j+1)-t_C(i,j) ;
        dx2= t_C(i,j+2)-t_C(i,j+1);
   
        dphidt(i,j) = (-(dx1/dx2)*phi(i,j+2) + ((dx1/dx2)+2+(dx2/dx1)) *phi(i,j+1) - (2+(dx2/dx1))*phi(i,j) )...
                                 / (dx1+dx2);
                             
%         dphidt(i,j) = (phi(i,j+1)-phi(i,j))/dx1;
        
        end
        
        if j==N-1                                                          % grad @ the Foot
        
        dx1= t_C(i,j)-t_C(i,j-1);
        dx2= t_C(i,j-1)-t_C(i,j-2);
   
        dphidt(i,j) = -(-(dx1/dx2)*phi(i,j-2) + ((dx1/dx2)+2+(dx2/dx1)) *phi(i,j-1) - (2+(dx2/dx1))*phi(i,j) )...
                                 / (dx1+dx2);
                             
%         dphidt(i,j) = (phi(i,j)-phi(i,j-1))/dx1;
                                    
        end
        
        if j~=1 && j~=N-1                                                  % grad in the centre
        
        dx1= t_C(i,j)-t_C(i,j-1);
        dx2= t_C(i,j+1)-t_C(i,j);
        
        dphidt(i,j)= ( phi(i,j+1)*(dx1/dx2) + ((dx2/dx1)-(dx1/dx2))*phi(i,j)  -   (dx2/dx1)*phi(i,j-1) ) ...
                                / (dx1+dx2);
                            
%         dphidt(i,j)= (phi (i,j+1) - phi (i,j-1) )/ (2*dx1);
                            
        end

    end
    
end
