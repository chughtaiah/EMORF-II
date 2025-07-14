function MM =ukf_ideal_self_modular_nsensors_xpp(y,x,x_0,dt,Q,R,ind)
% [y_out,xout,xout_0,dt_out,Q_out,R_out,ind_sima]=target_tracking(lambda_contam);

    m = x_0;
    P = Q;

    dim=size(y,1);
    n = size(m,1);
    alpha = 1;
    beta = 2;
    kappa = 0;

    lambda = alpha^2 * (n + kappa) - n;        

    [WM,WC]=weights(n,alpha,beta,lambda);


    MM = zeros(size(m,1), size(y,2) );
    PP = zeros(size(P,1),size(P,2),size(y,2) );
%     N=length(y)
  N=size(y,2);
    
    for k=1:N
        SX=sigma_gen(m,P,n,lambda);        
        % Propagate through the dynamic model 
        HX=dynamic_model(SX,dt);        
        
        % Compute the predicted mean and covariance
        [m,P]=mean_cov(WM,WC,HX);   
        P = P + Q;
        
        
        SX=sigma_gen(m,P,n,lambda);        
        HY=meas_model(SX,dim);   
        
        [mu,Uk]=mean_cov(WM,WC,HY);    
        C=cross_cov(SX,HY,WC,m,mu);
        
       
    
        % Compute the gain and updated mean and covariance  
        
        Rkdi=inver(R,ind(:,k));
        invS=Rkdi-Rkdi*((eye(dim)+Uk*Rkdi)^-1)*Uk*Rkdi;
        K=C*invS;
        

        m = m + K*(y(:,k) - mu);
        P = P-C*K';
    
        MM(:,k) = m;
        PP(:,:,k) = P;
        
        
        
    end  


%     time_mse_ukf=mean((x(1,:)-MM(1,:)).^2+(x(3,:)-MM(3,:)).^2);
%     mean(mean((x(3,:)-MM(3,:)).^2));
    time_mse_ukf=mean(mean((x(:,:)-MM(:,:)).^2));
%     plot(x(3,:))
%     hold on
%     plot(MM(3,:))
%     figure
%     plot(((x(3,:)-MM(3,:)).^2));
%     (trace(  ( (x(:,:)-MM(:,:) )'*( x(:,:)-MM(:,:) )  ) ))
%     a
    function SS=sigma_gen(x,P,n,lambda)  
    P=nearestSPD(P);
    A = chol(P,'lower');
    SS = [zeros(size(x)) A -A];
    SS = sqrt(n + lambda)*SS + repmat(x,1,size(SS,2));
end

function HX=dynamic_model(SX,dt)
     HX=zeros(size(SX,1),size(SX,2));
        for i=1:size(SX,2)
        w_t=SX(5,i);
        F_t=[1 sin(w_t*dt)/w_t 0 (cos(w_t*dt)-1)/w_t 0;...
        0 cos(w_t*dt) 0 -sin(w_t*dt) 0;...
        0 (1-cos(w_t*dt))/w_t 1 sin(w_t*dt)/w_t 0;...
        0 sin(w_t*dt) 0 cos(w_t*dt) 0;...
        0 0 0 0 1];
%         F_t=[1 1 0 0;...
%         0 1 0 0;...
%         0 0 1 1;...
%         0 0 0 1];
        HX(:,i) = F_t*SX(:,i);
        end
end

function [x,P]=mean_cov(WM,WC,HX)
%         x= zeros(size(xin));
%         P = zeros(size(Pin));
x=0;P=0;
        for i=1:size(HX,2)
            x = x + WM(i) * HX(:,i);
        end
        for i=1:size(HX,2)
            P = P + WC(i) * (HX(:,i) - x) * (HX(:,i) - x)';
        end
end

function [C]=cross_cov(SX,HY,WC,xm,mu)
    C  = zeros(size(SX,1),size(HY,1));
    for i=1:size(SX,2)
            C = C + WC(i) * (SX(:,i) - xm) * (HY(:,i) - mu)';
    end
end

% function HY=meas_model(SX,m)
%            for ii=1:size(SX,2)
%            for jj=1:m
%               HY(jj,ii)=sqrt((SX(1,ii)-(jj-1)*350)^2+(SX(3,ii)-(350*mod(jj+1,2)))^2);
%            end
%            end
%  
% end

function HY=meas_model(SX,m)
           for jj=1:m
           for ii=1:size(SX,2)
              HY(jj,ii)=sqrt((SX(1,ii)-(1-1)*350)^2+(SX(3,ii)-(350*mod(1+1,2)))^2)-sqrt((SX(1,ii)-(jj+1-1)*350)^2+(SX(3,ii)-(350*mod(jj+1+1,2)))^2) ;
           end
           end
 
end

function [WM,WC]=weights(n,alpha,beta,lambda)
    WM = zeros(2*n+1,1);
    WC = zeros(2*n+1,1);
    
    for j=1:2*n+1
        if j==1
            wm = lambda / (n + lambda);
            wc = lambda / (n + lambda) + (1 - alpha^2 + beta);
        else
            wm = 1 / (2 * (n + lambda));
            wc = wm;
        end
        WM(j) = wm;
        WC(j) = wc;
    end
end


function [inv_sig]=inver(Sig,It)
    inv_sig=zeros(size(Sig));
    indc_one=find(It==1);
    inv_sig(indc_one,indc_one)=Sig(indc_one,indc_one)^-1;
end
    
end
    