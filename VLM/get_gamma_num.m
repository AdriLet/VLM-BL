function [gamma_chord,gamma_span ] = get_gamma_num(GAMMA)
% 9-3-2011:  Rename the function to avoid conflict with the ITAT

M=length(GAMMA(:,1))+1;
N=length(GAMMA(1,:))+1;

gamma_chord=zeros(M-1,N);
gamma_span=zeros(M,N-1);

gamma_chord(:,1)=GAMMA(:,1);
gamma_chord(:,end)=-GAMMA(:,end);

for j=1:N-2   
    for i=1:M-1
        gamma_chord(i,j+1) = - (GAMMA(i,j) - GAMMA(i,j+1));   
    end   
end

gamma_span(1,:)= - GAMMA(1,:); 

for j=1:N-1
    for i=1:M-2
        gamma_span(i+1,j) = GAMMA(i,j) - GAMMA(i+1,j);    
    end 
end