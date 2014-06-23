function [X_V,Y_V,Z_V] = get_FcosSp (X,Y,Z,M,N)

M_inputs =length(X(:,1));
N_inputs =length(X(1,:));

%1 Origin at the tack of the sail

X0=X-X(1,1);
Y0=Y-Y(1,1);
Z0=Z-Z(1,1);

%2 Unheeled the datas (rotation around x)

theta = atan (Z0(1,end)/Y0(1,end));

Y0_rot =  Y0 * cos(theta) + Z0 * sin (theta);
Z0_rot =  - Y0 * sin(theta) + Z0 * cos (theta);


%3 Origin at the tack  + Adim (X and Y_rot) 

s = zeros(M_inputs,N_inputs);
t = zeros(M_inputs,N_inputs);

test=0;     %if pin head
if (X0(end,end)-X0(1,end))==0
    
    for i=1:M_inputs
        
        X0(i,end)=X0(1,end)+(i-1)*0.001/M_inputs;
    
    end
        
end


for j=1:N_inputs
    
    s (:,j) = (X0(:,j)-X0(1,j)) / (X0(end,j)-X0(1,j));
        
end

for i=1:M_inputs
    
    t (i,:) = (Y0_rot(i,:)-Y0_rot(i,1)) / (Y0_rot(i,end)-Y0_rot(i,1));                          

end


%4 Full cosine spacing in s,t coordinates (points on the Vortex Line (V)
%and the control Points (C)
 
for i=1:M
 
    for j=1:N

        s_V_FcosSp(i,j)= 0.5*(1-cos((i-1)*pi/(M-1)));  
        t_V_FcosSp(i,j)= 0.5*(1-cos((j-1)*pi/(N-1))); 
           
    end
    
end

    
%5 Interpolation in s and t coordinates

      
X_V_0rot = griddata(s,t,X0,s_V_FcosSp,t_V_FcosSp,'linear');
Y_V_0rot = griddata(s,t,Y0_rot,s_V_FcosSp,t_V_FcosSp,'linear');
Z_V_0rot = griddata(s,t,Z0_rot,s_V_FcosSp,t_V_FcosSp,'linear');   


%6 re-heel the datas

X_V_0 = X_V_0rot;
Y_V_0 = Y_V_0rot * cos (theta) - Z_V_0rot * sin (theta); 
Z_V_0 = Y_V_0rot * sin (theta) + Z_V_0rot * cos (theta); 
 

%7 Move the origin

X_V = X_V_0+X(1,1);
Y_V = Y_V_0+Y(1,1);
Z_V = Z_V_0+Z(1,1);


if test==1;     %if pin head
    
   X_V(:,end) = X(:,end);
   Y_V(:,end) = Y(:,end);
   Z_V(:,end) = Z(:,end);
    
end
