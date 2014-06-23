function [dCp_def_C] = get_deformed_dCp(sail,k0,k1,k2)

% Parameter to define a pressure distribution

load coef_middle.mat;

%% 1 Inputs

X_V=sail.X_V;
Y_V=sail.Y_V;
Z_V=sail.Z_V;

X_C=sail.X_C;
Y_C=sail.Y_C;
Z_C=sail.Z_C;

dCp_C=sail.dCp_C;

N=length(X_V(1,:));
M=length(X_V(:,1));

deg_pol=15;
Uinf=1;


%% 2 Find the glauert coeficient of the pressure distribution


% 2.1 Adim

% 2.1.1 Unheel the sail
    
phi = atan (Z_V (1,end)/Y_V(1,end));
        
Y_V_rot = Y_V * cos(phi) + Z_V * sin (phi);
Z_V_rot = -Y_V * sin(phi) + Z_V * cos (phi);

Y_C_rot = Y_C * cos(phi) + Z_C * sin (phi);
Z_C_rot = -Y_C * sin(phi) + Z_C * cos (phi);

for j=1:N-1
    
    % 2.1.1.2 Axis at the leading edge
    
    X_LE=[(X_V(1,j)+X_V(1,j+1))/2];
    Z_LE=[(Z_V_rot(1,j)+Z_V_rot(1,j+1))/2];
    
    X_TE=[(X_V(end,j)+X_V(end,j+1))/2];
    Z_TE=[(Z_V_rot(end,j)+Z_V_rot(end,j+1))/2];
    
    X_C_0(:,j)=X_C(:,j)-X_LE;
    Z_C_rot_0(:,j)=Z_C_rot(:,j)-Z_LE;

    % 2.1.1.3 At evey height the camber is untwisted
          
    twist(1,j)=atan((Z_TE-Z_LE) / (X_TE-X_LE));          
        
    X_C_0_rot (:,j)= X_C_0 (:,j)* cos(twist(1,j)) + Z_C_rot_0(:,j) * sin (twist(1,j));
    Z_C_rot_0_rot (:,j) = -X_C_0(:,j) * sin(twist(1,j)) + Z_C_rot_0(:,j) * cos (twist(1,j));
    
    % 2.1.1.4 [-]
    
    chord= ((Z_TE-Z_LE)^2 + (X_TE-X_LE)^2)^0.5;
    
    X_C_adim (:,j)= (X_C_0_rot (:,j)/ chord)-0.5;
    Z_C_adim (:,j)= Z_C_rot_0_rot (:,j)/ chord;
    
    teta (:,j) = acos (2*X_C_adim(:,j));
    
end



% 2.2 Fit the Glauert series using a least square approach

for j=1:N-1
                      
    %2.2.1 Definition of B
           
    for k=1:M-1
                
        B (k,1) = - (4/Uinf) * tan(teta (k,j)/2);
    
        for l=2:deg_pol
                                                
            B(k,l) = - (4/Uinf) * sin ((l-1) * teta (k,j));
                                                  
        end

    end

    A(:,j)=pinv(B)*dCp_C(:,j);
       

    %2.2.2 Relation between the Fourier and the Polynomial coeficients   

    for k=1:deg_pol                                                
    
        for l=1:deg_pol
                                        
            coef(k,l)=coef_middle(k,l);
                                
        end

    end

    %2.2.3 Constraint between the polynomial coeficients (to ensure Z(0) and
    %Z(1) = 0 

    w=zeros(deg_pol,1);
    for i=1:deg_pol                                                            % If odd, put 1
                                
        if (i/2)~=round(i/2)
                                       
            w(i,1)=1/(2^i);
                                   
        end

    end


    %2.2.4 Least square approach using the Moore-Penrose pseudoinverse Matrix
    %to find the polynomial coeficient of the camber (c)
           
    b=[coef w; w' 0];           
    d = [A(:,j) ; 0];           
    pseudo_inverse=pinv(b);                                                    %Moore-Penrose pseudoinverse Matrix (pinv).           
    e=pinv(b)*d;

    C_without_C0 = e(1:deg_pol);      

    
    %2.2.5 The constraint Glauert terms
    
    
    Gl(:,j)=coef*C_without_C0;

    
    clear C_without_C0
    clear coef
    clear e
    clear b
    clear d
    clear w
    

%% 3 Parameters 


Ysp=[Y_C_rot(1,1);Y_C_rot(1,round(N*0.5));Y_C_rot(1,end)];

facti(1,1:N-1)=polyval(polyfit(Ysp,[1;k0;1],2),Y_C_rot(1,:));
facti(2,:)=polyval(polyfit(Ysp,[1;k1;1],2),Y_C_rot(1,:));
facti(3,:)=polyval(polyfit(Ysp,[1;k2;1],2),Y_C_rot(1,:));
facti(4:deg_pol,:)=1;

% figure (33)
% hold on
% plot (facti(1,:),'b')
% plot (facti(2,:),'or')
% plot (facti(3,:),'*b')
% hold off

%% 4 Generate the deformed pressure distrubtion

% 4.1 dCp_def

   for k=1:M-1
                
       temp (1,j) = - (4/Uinf) * tan(teta (k,j)/2)*A(1,j)*facti(1,j);
  
       for l=2:deg_pol
                                                 
           temp(l,j) = - (4/Uinf) * sin ((l-1) * teta (k,j))*A(l,j)*facti(l,j);
                                                  
       end
       
       dCp_def_C(k,j)=sum(temp(:,j));

   end

end


