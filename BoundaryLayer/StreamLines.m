
function result=StreamLines(Sails)

           % figure
            %colormap jet
           % hold on
           % mesh( Main.X, Main.Y, Main.Z);
           % mesh( Head.X, Head.Y, Head.Z);
            %colorbar
           % colormap([0 0 0]);
           % axis 'equal';
%             axes=caxis;
%             caxis([axes(1,1) axes(1,2)]);
%             xlabel('X [m]')
%             ylabel('Y [m]')
%             zlabel('Z [m]')
%             hidden off

    for k=1:max(size(Sails));

                parfor j=2:(size(Sails(k).X_V,2)-1);
                    
%                     Sails(k).StreamLines.Ex=[ getStreamLine4(Sails(k),j,0) ; Sails(k).StreamLines.Ex ];
%                     Sails(k).StreamLines.In=[ getStreamLine4(Sails(k),j,1)  ; Sails(k).StreamLines.In ];
%                     Sails(k).StreamLines.Vmean=[ getStreamLine4(Sails(k),j,1/2)  ; Sails(k).StreamLines.Vmean ];

                            Up(j-1,1)=getStreamLine4(Sails(k),j,0) ;
                           Lo(j-1,1)=getStreamLine4(Sails(k),j,1) ;
                            Vmean(j-1,1)=getStreamLine4(Sails(k),j , 1/2) ;
                end
                
                
                     Sails(k).StreamLines.Up=Up;
                    Sails(k).StreamLines.Lo=Lo;
                    Sails(k).StreamLines.Vmean=Vmean;
                    
                    clearvars Ex  In  Vmean ; 

    end
         result=Sails;
   end
