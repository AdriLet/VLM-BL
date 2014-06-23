function [B] = get_InflCoefs (sailInfling,sailInfled)
    %matlabpool;
    X_V=sailInfling.X_V;
    Y_V=sailInfling.Y_V;
    Z_V=sailInfling.Z_V;

    wake_X= sailInfling.wake_X;
    wake_Y= sailInfling.wake_Y;
    wake_Z= sailInfling.wake_Z;

    X_C=sailInfled.X_C;
    Y_C=sailInfled.Y_C;
    Z_C=sailInfled.Z_C;

    nx= sailInfled.nx;
    ny= sailInfled.ny;
    nz= sailInfled.nz;

    M_infled=length (X_C(:,1))+1;
    N_infled=length (X_C(1,:))+1;

    M_infling=length (sailInfling.X_V(:,1));
    N_infling=length (sailInfling.X_V(1,:));

    B=zeros((M_infled-1)*(N_infled-1),(M_infling-1)*(N_infling-1));
    temp=1;

    b_panel_u=zeros(M_infling-1,N_infling-1);
    b_panel_v=zeros(M_infling-1,N_infling-1);
    b_panel_w=zeros(M_infling-1,N_infling-1);

    for l=1:N_infled-1
        for k=1:M_infled-1

            for j=1:N_infling-1
               for i=1:M_infling-2;
                
                        [b1_u,b1_v,b1_w] = get_V_from_vortl (X_C(k,l),Y_C(k,l),Z_C(k,l),     X_V(i,j), Y_V(i,j),Z_V(i,j),                 X_V(i+1,j),Y_V(i+1,j),Z_V(i+1,j),        1);
                        [b2_u,b2_v,b2_w] = get_V_from_vortl (X_C(k,l),Y_C(k,l),Z_C(k,l),     X_V(i+1,j),Y_V(i+1,j),Z_V(i+1,j),           X_V(i+1,j+1),Y_V(i+1,j+1),Z_V(i+1,j+1),  1);
                        [b3_u,b3_v,b3_w] = get_V_from_vortl (X_C(k,l),Y_C(k,l),Z_C(k,l),     X_V(i+1,j+1),Y_V(i+1,j+1),Z_V(i+1,j+1),     X_V(i,j+1),Y_V(i,j+1),Z_V(i,j+1),        1);
                        [b4_u,b4_v,b4_w] = get_V_from_vortl (X_C(k,l),Y_C(k,l),Z_C(k,l),     X_V(i,j+1),Y_V(i,j+1),Z_V(i,j+1),           X_V(i,j),Y_V(i,j),Z_V(i,j),              1);

                        b_panel_u(i,j)=b1_u + b2_u + b3_u + b4_u;
                        b_panel_v(i,j)=b1_v + b2_v + b3_v + b4_v;
                        b_panel_w(i,j)=b1_w + b2_w + b3_w + b4_w;

                        %    4
                        %    _____    ->Y
                        %   |          |
                        %1 |          |3
                        %   |_____|
                        %    2
                        %\|/ X
                   
                end

                [b1_u,b1_v,b1_w] = get_V_from_vortl (X_C(k,l),Y_C(k,l),Z_C(k,l),     X_V(end-1,j),Y_V(end-1,j),Z_V(end-1,j),               X_V(end,j),Y_V(end,j),Z_V(end,j),                  1);
                [b3_u,b3_v,b3_w] = get_V_from_vortl (X_C(k,l),Y_C(k,l),Z_C(k,l),     X_V(end,j+1),Y_V(end,j+1),Z_V(end,j+1),               X_V(end-1,j+1),Y_V(end-1,j+1),Z_V(end-1,j+1),      1);
                [b4_u,b4_v,b4_w] = get_V_from_vortl (X_C(k,l),Y_C(k,l),Z_C(k,l),     X_V(end-1,j+1),Y_V(end-1,j+1),Z_V(end-1,j+1),         X_V(end-1,j),Y_V(end-1,j),Z_V(end-1,j),            1);

                [b1inf_u, b1inf_v, b1inf_w] = get_V_from_semiinfvortl (X_C(k,l),Y_C(k,l),Z_C(k,l),       X_V(end,j),Y_V(end,j),Z_V(end,j),         wake_X(2,j),wake_Y(2,j),wake_Z(2,j),        1);
                [b3inf_u, b3inf_v, b3inf_w] = get_V_from_semiinfvortl (X_C(k,l),Y_C(k,l),Z_C(k,l),       X_V(end,j+1),Y_V(end,j+1),Z_V(end,j+1),   wake_X(2,j+1),wake_Y(2,j+1),wake_Z(2,j+1),      1);

                b_panel_u(M_infling-1,j)=b1_u + b3_u + b4_u + b1inf_u - b3inf_u;
                b_panel_v(M_infling-1,j)=b1_v + b3_v + b4_v + b1inf_v - b3inf_v;
                b_panel_w(M_infling-1,j)=b1_w + b3_w + b4_w + b1inf_w - b3inf_w;

            end

            [b_panel] = get_projection (b_panel_u,b_panel_v,b_panel_w,nx(k,l),ny(k,l),nz(k,l));

            B(temp,:)=reshape (b_panel,1,(M_infling-1)*(N_infling-1));
            temp=temp+1;

        end
        
    end
   % parpool close;
end