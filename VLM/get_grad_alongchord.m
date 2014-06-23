function [dphids]= get_grad_alongchord(s_C,t_C,phi)

M=length(s_C(:,1))+1;
N=length(s_C(1,:))+1;

% dphi/ds
%--------

for j=1:N-1
    
    for i=1:M-1
        
        
        if i==1                                                            % grad @ panel of the LE
            
        dx1=s_C(i+1,j)-s_C(i,j);
        dx2=s_C(i+2,j)-s_C(i+1,j);
   
        dphids(i,j) = (-(dx1/dx2)*phi(i+2,j) + ((dx1/dx2)+2+(dx2/dx1)) *phi(i+1,j) - (2+(dx2/dx1))*phi(i,j) )...
                                 / (dx1+dx2);

%         dphids(i,j) = (phi(i+1,j)-phi(i,j))/dx1;
        
        end
        
        if i==M-1                                                          % grad @ panel of the TE
        
        dx1= s_C(i,j)-s_C(i-1,j);
        dx2= s_C(i-1,j)-s_C(i-2,j);
   
        dphids(i,j) = -(-(dx1/dx2)*phi(i-2,j) + ((dx1/dx2)+2+(dx2/dx1)) *phi(i-1,j) - (2+(dx2/dx1))*phi(i,j) )...
                                 / (dx1+dx2);

%         dphids(i,j) = (phi(i,j)-phi(i-1,j))/dx1;
                                    
        end
        
        if i~=1 && i~=M-1                                                  % grad in the centre
        
        dx1= s_C(i,j)-s_C(i-1,j);
        dx2= s_C(i+1,j)-s_C(i,j);
        
        dphids(i,j)= ( phi(i+1,j)*(dx1/dx2) + ((dx2/dx1)-(dx1/dx2))*phi(i,j)  -   (dx2/dx1)*phi(i-1,j) ) ...
                                / (dx1+dx2);

%         dphids(i,j)= (phi (i+1,j) - phi (i-1,j) )/ (2*dx1);
                            
      
        end


    end
    
end
