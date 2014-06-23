function [s,t] = get_sandt (X,Y,Z)

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
