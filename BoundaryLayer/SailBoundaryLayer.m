function [ Out_Sail] = SailBoundaryLayer( Sail)
%SAILBOUNDARYLAYER  This function calculates the boundary layer on a whole sail. Skin friction and stall position are calculated (2D) along streamlines.



%%  Extrado
 
M=size(Sail.StreamLines.Up,1);
    for k=[1:M]

        
        X=Sail.StreamLines.Up(k).X;
        Y=Sail.StreamLines.Up(k).Y;
        Z=Sail.StreamLines.Up(k).Z;
        V=Sail.StreamLines.Up(k).V;
        Vinf=Sail.StreamLines.Up(k).Vinf;
        N=max(size(V));
        dS(2:N)=sqrt((X(1:N-1)-X(2:N)).^2+(Y(1:N-1)-Y(2:N)).^2+(Z(1:N-1)-Z(2:N)).^2);          % S is curvilinear coordinate.
        dS(1)=0;
        S(1)=0;
        for i=2:N
            S(i)=S(i-1)+dS(i);
        end
        
        if (size(S)==size(V))
        else
            disp('ho')
        end 
        
        
        
      %  [ detach_lam, transition, detach_turb, X_tr, X_det, Cf] =BL_detachment_old(S,V,Vinf);
        [ detach_turb, Xstall, Ystall, Zstall,Cf,H,theta] =HeadsMethod(Sail.StreamLines.Up(k,1));
        if detach_turb
            Sail.StreamLines.DetachLine.Up.X(k)=Xstall;
            Sail.StreamLines.DetachLine.Up.Y(k)=Ystall;
            Sail.StreamLines.DetachLine.Up.Z(k)=Zstall;
        else
            Sail.StreamLines.DetachLine.Up.X(k)=NaN;
            Sail.StreamLines.DetachLine.Up.Y(k)=NaN;
            Sail.StreamLines.DetachLine.Up.Z(k)=NaN;
        end
                Sail.StreamLines.Up(k).H=H;
                Sail.StreamLines.Up(k).theta=theta;
                Sail.StreamLines.Up(k).Cf=Cf;
                clearvars S X Y Z theta H Cf;
    end
    
%% Intrado
M=size(Sail.StreamLines.Lo,1);
   for k=[1:M]
        
        
        X=Sail.StreamLines.Lo(k).X;
        Y=Sail.StreamLines.Lo(k).Y;
        Z=Sail.StreamLines.Lo(k).Z;
        V=Sail.StreamLines.Lo(k).V;
        Vinf=Sail.StreamLines.Lo(k).Vinf;
        N=max(size(V));
        dS(2:N)=sqrt((X(1:N-1)-X(2:N)).^2+(Y(1:N-1)-Y(2:N)).^2+(Z(1:N-1)-Z(2:N)).^2);          % S is curvilinear coordinate.
        dS(1)=0;
        S(1)=0;
        for i=2:N
            S(i)=S(i-1)+dS(i);
        end
      %  [ detach_lam, transition, detach_turb, X_tr, X_det,Cf] =BL_detachment_old(S,V,Vinf);
          [ detach_turb, Xstall, Ystall, Zstall,Cf,H,theta] =HeadsMethod(Sail.StreamLines.Lo(k,1));
        if detach_turb
            Sail.StreamLines.DetachLine.Lo.X(k)=Xstall;
            Sail.StreamLines.DetachLine.Lo.Y(k)=Ystall;
            Sail.StreamLines.DetachLine.Lo.Z(k)=Zstall;
        else
            Sail.StreamLines.DetachLine.Lo.X(k)=NaN;
            Sail.StreamLines.DetachLine.Lo.Y(k)=NaN;
            Sail.StreamLines.DetachLine.Lo.Z(k)=NaN;   
            
        end
                Sail.StreamLines.Lo(k).H=H;
                Sail.StreamLines.Lo(k).theta=theta;
                Sail.StreamLines.Lo(k).Cf=Cf;
                clearvars S X Y Z theta H Cf;
    end
    


Out_Sail=Sail;

    
end

