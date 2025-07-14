function [y,x,x_0,dt] = target_tracking_var_dim_nsensors(dim,R,Q,ind_sim,T,alpha_contam,sig)
% close all;
dt=1; % Sampling Time

% Process Noise Covariance


% State at t_0
%   x_0 = [-10000;10;5000;-5;-0.053];
% x_0 = mvnrnd([-10000;10;5000;-5;-3*pi/180],Q)';
x_0 = mvnrnd([0;1;0;-1;-3*pi/180],Q)';



w_t=x_0(5);


x_t=x_0;
%probabilities of occurence of outlier in a given dimension; %level of contamination 

% Running Simulation Till time T




for t=1:T
    % Matrix for Function F at time T  
    F_t=[1 sin(w_t*dt)/w_t 0 (cos(w_t*dt)-1)/w_t 0;...
        0 cos(w_t*dt) 0 -sin(w_t*dt) 0;...
        0 (1-cos(w_t*dt))/w_t 1 sin(w_t*dt)/w_t 0;...
        0 sin(w_t*dt) 0 cos(w_t*dt) 0;...
        0 0 0 0 1]; 

%     F_t=[1 1 0 0;...
%         0 1 0 0;...
%         0 0 1 1;...
%         0 0 0 1]; 
    
    % Generating Zero Mean Gaussian Process Noise with Covariance Q
    
    ut = mvnrnd([0 0 0 0 0],Q);  
%       ut = mvnrnd([0 0 0 0],Q);  
  
    % Obtaining State at time t
    x(:,t)=F_t*x_t+ut';
    
    x_t=x(:,t);
    w_t=x_t(5); 
     
    
%     samples_R=mvnrnd(zeros(dim,1),R);

        
    for ii=1:dim
        TOA(ii,t)=sqrt((x_t(1)-(ii-1)*350)^2+(x_t(3)-(350*mod(ii+1,2)))^2)+mvnrnd(0,sig(ii));
    end

            
    for ii=1:dim-1
        y(ii,t)=TOA(1,t)-TOA(ii+1,t);
    end
            
    for ii=1:dim-1
            if(ind_sim(ii,t)==0)
            y(ii,t)=y(ii,t)+mvnrnd(0,R(ii,ii)*alpha_contam);
            end
    end


end


% function Ui = randU(n)
% Xi = (randn(n) )/sqrt(2);
% [Qi,Ri] = qr(Xi);
% Ri = diag(diag(Ri)./abs(diag(Ri)));
% Ui = Qi*Ri;
% end


end
