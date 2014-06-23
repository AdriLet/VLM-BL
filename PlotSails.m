function PlotSails(Sails)


figure(1)
clf
%     subplot(1,2,1);
    colormap jet
    hold on
    for i=1:max(size(Sails))
        mesh(Sails(i).X, Sails(i).Y, Sails(i).Z);
    end

    %colorbar
    colormap([0 0 0]);
    axis 'equal';
    axes=caxis;
    caxis([axes(1,1) axes(1,2)]);
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    hidden off

    
    for k=1:max(size(Sails))
        
        hold on
        for i=1:max(size(Sails(k).StreamLines.Up))
            line=Sails(k).StreamLines.Up;
            plot3(line(i,1).X,line(i,1).Y,line(i,1).Z,'-r','linewidth',1);
            line=Sails(k).StreamLines.Lo;
            plot3(line(i,1).X,line(i,1).Y,line(i,1).Z,'-g','linewidth',1);
%             line=Sails(k).StreamLines.Vmean;
%             plot3(line(i,1).X,line(i,1).Y,line(i,1).Z,'-b','linewidth',2);
        end
        
        
        plot3(Sails(k).StreamLines.DetachLine.Up.X,Sails(k).StreamLines.DetachLine.Up.Y,Sails(k).StreamLines.DetachLine.Up.Z,'-+r','linewidth',3);
        plot3(Sails(k).StreamLines.DetachLine.Lo.X,Sails(k).StreamLines.DetachLine.Lo.Y,Sails(k).StreamLines.DetachLine.Lo.Z,'-+g','linewidth',3);
        
    end
    
  title(strcat('AWA =  ',num2str(Sails(1).AWA), '^{\circ} ;  TWS = ', num2str(Sails(1).TWS) ,   ' m/s ; BoatSpeed =  ', num2str(Sails(1).Bs), ' m/s' ) ,'FontWeight','bold')
hold off


print(figure(1),'-dpdf',strcat('AWA=',num2str(Sails(1).AWA),'.pdf'))
clf
% subplot(1,2,2);
% % 6.2 dCp
% %--------
% 
% % 6.1 Shapes + wakes
% %--------------------hold on
% hold on
%     for k=1:max(size(Sails))
% 
% colormap
% mesh(Sails(k).X,Sails(k).Y,Sails(k).Z,Sails(k).dCp_C);
% plot3(Sails(k).wake_X,Sails(k).wake_Y,Sails(k).wake_Z,'-r');
% %colormap jet
% %surf(Sails(k).X_C,Sails(k).Y_C,Sails(k).Z_C,Sails(k).dCp_C);
% view(0,90)
% axis('equal')
% xlabel('X [m]')
% ylabel('Y [m]')
% zlabel('Z [m]')
%     end
%     hold off
end 


