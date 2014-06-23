function [Vx_ind,Vy_ind,Vz_ind ] = get_Vind(sailInfling,sailInfled)

X_V=sailInfling.X_V;
Y_V=sailInfling.Y_V;
Z_V=sailInfling.Z_V;

wake_X= sailInfling.wake_X;
wake_Y= sailInfling.wake_Y;
wake_Z= sailInfling.wake_Z;

X_C=sailInfled.X_C;
Y_C=sailInfled.Y_C;
Z_C=sailInfled.Z_C;


M_infled=length (sailInfled.X_C(:,1))+1; 
N_infled=length (sailInfled.X_C(1,:))+1; 

M_infling=length (sailInfling.X_V(:,1));
N_infling=length (sailInfling.X_V(1,:));

% 1 The gamma
[gamma_chord,gamma_span ] = get_gamma_num(sailInfling.GAMMA);

% 2 Computation of the induced velocity.
% temp=1;

Vx_ind_fromchordVL=zeros(M_infling-1,N_infling);
Vy_ind_fromchordVL=zeros(M_infling-1,N_infling);
Vz_ind_fromchordVL=zeros(M_infling-1,N_infling);

Vx_ind_fromspanVL=zeros(M_infling,N_infling-1);
Vy_ind_fromspanVL=zeros(M_infling,N_infling-1);
Vz_ind_fromspanVL=zeros(M_infling,N_infling-1);

Vx_ind_fromwake=zeros(1,N_infling);
Vy_ind_fromwake=zeros(1,N_infling);
Vz_ind_fromwake=zeros(1,N_infling);

Vx_ind=zeros(M_infled-1,N_infled-1);
Vy_ind=zeros(M_infled-1,N_infled-1);
Vz_ind=zeros(M_infled-1,N_infled-1);

for l=1:N_infled-1
    for k=1:M_infled-1
              
        for j=1:N_infling
            for i=1:M_infling-1 
                [Vx_ind_fromchordVL(i,j),Vy_ind_fromchordVL(i,j),Vz_ind_fromchordVL(i,j)] = get_V_from_vortl (X_C(k,l),Y_C(k,l),Z_C(k,l),            X_V(i,j),Y_V(i,j),Z_V(i,j),         X_V(i+1,j),Y_V(i+1,j),Z_V(i+1,j),gamma_chord(i,j));     
            end       
        end
        
        for j=1:N_infling-1
            for i=1:M_infling 
                [Vx_ind_fromspanVL(i,j),Vy_ind_fromspanVL(i,j),Vz_ind_fromspanVL(i,j)] = get_V_from_vortl (X_C(k,l),Y_C(k,l),Z_C(k,l),          X_V(i,j),Y_V(i,j),Z_V(i,j),         X_V(i,j+1),Y_V(i,j+1),Z_V(i,j+1),gamma_span(i,j));
            end
        end
         
       for j=1:N_infling  
         [Vx_ind_fromwake(j),Vy_ind_fromwake(j),Vz_ind_fromwake(j)] = get_V_from_semiinfvortl (X_C(k,l),Y_C(k,l),Z_C(k,l),           X_V(end,j),Y_V(end,j),Z_V(end,j),         wake_X(2,j),wake_Y(2,j),wake_Z(2,j),gamma_chord(end,j));  
       end

       Vx_ind(k,l) = sum (sum (Vx_ind_fromchordVL)) + sum (sum (Vx_ind_fromspanVL))+ sum(Vx_ind_fromwake);
       Vy_ind(k,l) = sum (sum (Vy_ind_fromchordVL)) + sum (sum (Vy_ind_fromspanVL))+ sum (Vy_ind_fromwake);
       Vz_ind(k,l) = sum (sum (Vz_ind_fromchordVL)) + sum (sum (Vz_ind_fromspanVL)) + sum (Vz_ind_fromwake);
  
    end
end
