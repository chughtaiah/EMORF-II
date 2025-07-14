function [R_out,Q_out,ind_sima,sig]=gen(dim,lambda_contam,T)

dt=1;

ind_sima=ones(dim-1,T);

% U=randU(dim);

sig=zeros(dim,1);

for ii=1:dim
sig(ii)=10;
end

R_out=sig(1)*ones(dim-1,dim-1);

for jj=1:dim-1
R_out(jj,jj)=R_out(jj,jj)+sig(jj+1);
end


q1=.1;
q2=1.75*(10^-4);


Q_out= [q1*dt^3/3 q1*dt^2/2 0 0 0;...
    q1*dt^2/2 q1*dt 0 0 0;...
    0 0 q1*dt^3/3 q1*dt^2/2 0;...
    0 0 q1*dt^2/2 q1*dt 0;...
    0 0 0 0 q2]; % This Matrix is defined only once because it is constant for all
%                  time points
                 
% Q_out= [q1*dt^3/3 q1*dt^2/2 0 0 ;...
%     q1*dt^2/2 q1*dt 0 0 ;...
%     0 0 q1*dt^3/3 q1*dt^2/2;...
%     0 0 q1*dt^2/2 q1*dt]; 

% Q_out=.1*eye(4);

% Q_out=[.5*dt^2 0;1*dt 0;0 .5*dt^2;0 1*dt]*.1*eye(2)*[.5*dt^2 0;1*dt 0;0 .5*dt^2;0 1*dt]';
                 
for t=1:T
    for ii=1:dim-1
    if rand<1-(1-lambda_contam)^2
    ind_sima(ii,t)=0;
    end
    end
end







end