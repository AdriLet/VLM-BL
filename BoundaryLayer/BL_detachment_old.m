function [ detach_lam, transition, detach_turb, X_tr, X_det,Cf] = BL_detachment_old( x, Ue, Uinf)
%Constantes du problème
%clear;
%clc;
%close all
laminar=false; %si ce paramètre vaut 1 alors la couche limite est liminaire sinon turbulent
laminar_method=0;        %choix du critère en laminaire. 0=Thwaites, 1=Stratford 


rho=1.184;% densité de l'air à 25° a 1 atm
nu=1.8e-5*rho;              %viscosité cinématique de l'air


N=size(x,2); %nombre de points de calculs selon alpha (nombre de points de calculs -1)
p=-0.5*rho*Ue.^2;               
Inte_Ue5=0;
Inte_Ue6=0;
Inte_Ue3=0;
pm=min(p);
m=min(find(p==pm));
xm=x(m);
Um=Ue(m);
Umax=max(Ue);
X_tr=0;
 
separated=false; %si ce paramètre vaut true alors la couche limite est séparée 
%Cp_=[(p(1)-pm)/(0.5*rho*Um^2)];


        
Cp_=ones(1,N)-Ue.^2/Um^2;
Cf=[];
Lambda=zeros(1,N);
strat=zeros(1,N);
%alpha=180/pi*acos(ones(1,size(x,2))-x*1/R);
detach_lam=false;
detach_turb=false;
transition=false;
N=size(x,2);
X_det=x(N);
sigma=13*ones(1,N);
Mitchel=zeros(1,N);



U_tr=1;
%Calcul de xm_ pour le critère de Startford en laminaire.
xm_=0;
for i=1:m-1
    
   xm_=xm_ +  Ue(i)^5*(x(i+1)-x(i));
end
xm_=xm_/Um^5;




%debut des itérations le long de la ligne de courant
for i=2:N-1
    
    
    if and(laminar,not(separated))
          %Cp_=[Cp_ (p(i)-pm)/(0.5*rho*Um^2)];
          
          Inte_Ue5=Inte_Ue5+Ue(i)^5*(x(i)-x(i-1));

      
      
            theta=sqrt(0.45*nu/Ue(i)^6*Inte_Ue5);
%Critère de Twaites en couche limite laminaire
      if laminar_method==0
              
              lambda=theta^2/nu*(Ue(i)-Ue(i-1))/(x(i)-x(i-1));
                Lambda(i)= lambda;
              if and(lambda<-0.095,i>m)
%                   disp('laminar separation according Thwaites at')
%                  alpha(i)    
                 separated=true;
                detach_lam=true;
                X_det=x(i);
               
              end
              
      end
%Calcul du coefficient de frottement selon Thwaites.
%             l=0;
%             if (0<lambda)
%                 l=0.22+1.57*lambda-1.8*lambda^2;
%             elseif (lambda<=0)
%                 l=0.22+1.402*lambda+0.018*lambda/(lambda+0.107);
%             end
%             Cf=[Cf l*nu/(0.5*Ue(i)*theta)];
%                     
%     

      
%Critère de Stratford en couche limite laminaire
        
            if and(laminar_method==1,x(i)>xm) 
                x_=x(i)-xm+xm_;
                dCp_=(Cp_(i)-Cp_(i-1))/(x(i)-x(i-1));
                Strat=Cp_(i)*(x_*dCp_)^2;
                strat(i)=Strat;
                if Strat>0.0104                  
%                      disp('laminar separation according Stratford at');
%                          alpha(i)  
                        separated=true;
                        detach_lam=true;
                         X_det=x(i);
                end
            end
      

%critère de Michell pour la transition laminaire-turbulent
    
          Re_theta=theta*Uinf/nu;
          Re_x=Uinf*x(i)/nu;
          Mitchel(i)=Re_theta / 1.174*(1+22.400/Re_x) *Re_x^(0.46);
          if Mitchel(i)<1
              laminar = false;                                    
%               disp('transition according Mitchel at ');
%               x(i)
              transition = true;
              X_tr=x(i);
              U_tr=Ue(i);
              %alpha(i)
              if X_tr>xm
                disp('Be carefull, Xm is before Xtr');
              end
          end
        end 
          


%critère de Straford en couche limite turbulente
    if and( not(laminar), not(separated))

             Re_x=Uinf*x(i)/nu;
             Cf=[Cf 0.0594/(Re_x^0.2)];
             %Cp_=[Cp_ (p(i)-pm)/(0.5*rho*Um^2)];
             dCp_=(Cp_(i)-Cp_(i-1))/(x(i)-x(i-1));

             if x(i)<xm
                Inte_Ue3=Inte_Ue3+Ue(i)^3*(x(i)-x(i-1));
             end
             
             if x(i)>xm
                 Xm_= 38.5*nu^(5/3)*U_tr^(-25/8)*Inte_Ue5^(5/8)+Um^(-3)*Inte_Ue3;
                 x_=x(i)-xm+Xm_;
                 Re_xm_=Uinf*Xm_/nu;
                 Strat=Cp_(i)*(x_*dCp_)^(0.5)/(Re_xm_/10^6)^(0.1);
                strat(i)=Strat;
                d2p=(p(i+1)-2*p(i)+p(i-1));    %on ne s'interresse qu'au signe, donc ne ne prend pas en compte le pas.
                 if d2p>0
                     if Strat>0.39
%                           disp('turbulent separation according Stratford at ')
%                           x(i)
                          separated=true;
                           detach_turb=true;
                           X_det=x(i);
                           if Cp_(i)>(4/7)
                            disp('Be careful Cp_ > 4/7 so criterium out of range of use');
                                 Cp_(i)
                             end
                     end
                 else
                     if Strat>0.35
%                           disp('turbulent separation according Stratford at ')
%                           x(i)
                          separated=true;
                          detach_turb=true;
                          X_det=x(i);
                          if Cp_(i)>(4/7)
                            disp('Be careful Cp_ > 4/7 so criterium out of range of use');
                            Cp_(i)
                          end
                     end
                 end
             
             end
       
             
             
    end  

 
end
       DUe=Ue(2:end)-Ue(1:end-1);
       DX=x(2:end)-x(1:end-1);
       Ue_prime=DUe./DX;
       DUe_prime=Ue_prime(2:end)-Ue_prime(1:end-1);
       Ue_prime_prime=DUe_prime./DX(1:end-1);
       sigma_bis=Ue(1:end-2).*Ue_prime_prime./(Ue_prime(1:end-1).^2);
% figure 
% plot(x/x(N),Cp_)
% hold off
% title('Cp(x)');
% xlabel('x/x_{max}');
% ylabel('Cp');
% 
% figure
% plot(x/x(N),strat)
% hold off
% title('Startford Criterion');
% xlabel('x/x_{max}');
% ylabel('Startford number');
% 
% 
% figure
% plot(x/x(N),Lambda)
% hold off
% 
% title('Twaites Method');
% xlabel('x/x_{max}');
% ylabel('lambda');
% 
% figure
%  plot(x/x(N),Mitchel)
%  hold off
% 
% title('Mitchel criterion');
% xlabel('x/x_{max}');
% ylabel('Mitchel number');
% 
% 
% figure
%  plot(x(1:end-2)/x(N),sigma_bis)
%  hold off
% 
% title('sigma');
% xlabel('x/x_{max}');
% ylabel('sigmabis');


end

