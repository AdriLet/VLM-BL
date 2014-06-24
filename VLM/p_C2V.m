function [p_V]=p_C2V(Sail,p_C)
X_V=Sail.X_V;
Y_V=Sail.Y_V;
Z_V=Sail.Z_V;
X_C=Sail.X_C;
Y_C=Sail.Y_C;
Z_C=Sail.Z_C;
M=length(X_V(:,1));
N=length(X_V(1,:));

% 1 Points at the middle the VL in spanwise direction
for i=1:M
    for j=1:N-1
        X_midSpanVL(i,j)= X_V(i,j)+ [X_V(i,j+1)- X_V(i,j)]/2;
        Y_midSpanVL(i,j)= Y_V(i,j)+ [Y_V(i,j+1)- Y_V(i,j)]/2;
        Z_midSpanVL(i,j)= Z_V(i,j)+ [Z_V(i,j+1)- Z_V(i,j)]/2;
    end
end

% 2 Control points to the middle of the spanwise vortex lines
% 2.1 Distance (along the chord) between the CP + Distance (along the chord) between the middle the VL .
for j=1:N-1
    c_CP2midSpanVL(:,j)= ((X_C(:,j)-X_midSpanVL(1,j)).^2 + (Y_C(:,j)-Y_midSpanVL(1,j)).^2 + (Z_C(:,j)-Z_midSpanVL(1,j)).^2).^0.5 ;
    c_midSpanVL2midSpanVL(:,j)= ((X_midSpanVL(:,j)-X_midSpanVL(1,j)).^2 + (Y_midSpanVL(:,j)-Y_midSpanVL(1,j)).^2 + (Z_midSpanVL(:,j)-Z_midSpanVL(1,j)).^2).^0.5 ;
end

% 2.2 Extrapolation
for j=1:N-1
    p_atmidSpanVL(:,j) = interp1(c_CP2midSpanVL(:,j),p_C(:,j),c_midSpanVL2midSpanVL(:,j),'cubic');
    
end

% 3 Middle of the spanwise vortex lines to the vortex lines
% 3.1 Distance (along the span) between the middle of the vortex lines + Distance (along the span)
%between the extremities of the VL
for i=1:M
    c_midSpanVL2VL(i,:)= ((X_midSpanVL(i,:)-X_V(i,1)).^2 + (Y_midSpanVL(i,:)-Y_V(i,1)).^2 + (Z_midSpanVL(i,:)-Z_V(i,1)).^2).^0.5 ;
    c_V(i,:)= ((X_V(i,:)-X_V(i,1)).^2 + (Y_V(i,:)-Y_V(i,1)).^2+ (Z_V(i,:)-Z_V(i,1)).^2).^0.5 ;
end

% 3.2 Extrapolation
for i=1:M
    p_V(i,:) = interp1(c_midSpanVL2VL(i,:),p_atmidSpanVL(i,:),c_V(i,:),'cubic');
end
end
