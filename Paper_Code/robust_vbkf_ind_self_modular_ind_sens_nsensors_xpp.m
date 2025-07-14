function xpp =robust_vbkf_ind_self_modular_ind_sens_nsensors_xpp(y,x,x_0,dt,Q,R)


N=size(y,2);
xe=x_0;
Pe=Q;

    n = size(xe,1);
    alpha = 1;
    beta = 2;
    kappa = 0;
    lambda = alpha^2 * (n + kappa) - n;
    
    [WM,WC]=weights(n,alpha,beta,lambda);
    
    xpp=[];
        e_0=0.9;
        f_0=0.1;
    dim=size(y,1);

for k=1:N
        d=1;
        %generate sigmas
        SX=sigma_gen(xe,Pe,n,lambda);        
        % Propagate through the dynamic model 
        HX=dynamic_model(SX,dt);        
        
        % Compute the predicted mean and covariance
        [xm,Pm]=mean_cov(WM,WC,HX);   
        Pm = Pm + Q;

        z_t=ones(dim,1);
        e_t = e_0*ones(dim,1);
        f_t = f_0*ones(dim,1);
    
        j=1;
    
        % Form sigma points for measurement step and
        % propagate throught the measurement model
        SX=sigma_gen(xm,Pm,n,lambda);        
        HY=meas_model(SX,dim);   
    
        
        [mu,Uk]=mean_cov(WM,WC,HY);    
        C=cross_cov(SX,HY,WC,xm,mu);
        
    Sig_zt=R;
    
    K = C*((Uk + Sig_zt)^-1);  

    
    xp(:,j)=xm+K*(y(:,k)-mu);
    Pp=Pm-C*K';

    SX=sigma_gen(xp(:,j),Pp,n,lambda);
    HY=meas_model(SX,dim);
        
    [my,Yk]=mean_cov(WM,WC,HY);              
    Bk=(my-y(:,k))*(my-y(:,k))'+Yk;
    
    
    for iz=1:dim
    
    
    z_t(iz)=1;
    pz1=exp( -.5*Bk(iz,iz)/R(iz,iz)+( psi(e_t(iz))-psi(e_t(iz)+f_t(iz))));
    z_t(iz)=0;
    pz0=exp( ( psi(f_t(iz))-psi(e_t(iz)+f_t(iz))) );
    rat=pz1/(pz1+pz0);
    
        if(rat>.5)
        z_t(iz)=1;
        else
        z_t(iz)=0;
        end
        
        if((z_t(iz))==0)
        Sig_zt(iz,iz)=inf;
        else
        Sig_zt(iz,iz)=R(iz,iz);
        end
    end 

    
    
    %% Updating e_t and f_t
    e_t = e_0 + z_t;
    f_t = f_0 + 1 - z_t; 

    
    while(d>10^-4)
    j=j+1;

    Rkdi=inver_sig_test(Sig_zt);
    invS=Rkdi-Rkdi*((eye(dim)+Uk*Rkdi)^-1)*Uk*Rkdi;
    K=C*invS;


    xp(:,j)=xm+K*(y(:,k)-mu);
    Pp=Pm-C*K';    

    SX=sigma_gen(xp(:,j),Pp,n,lambda);
    HY=meas_model(SX,dim);        
    [my,Yk]=mean_cov(WM,WC,HY);  
        
    
    Bk=(my-y(:,k))*(my-y(:,k))'+Yk;
%     ln_pi_expected = psi(e_t)-psi(e_t+f_t);
%     ln_one_minus_pi_expected=psi(f_t)-psi(e_t+f_t);
%     z_t = (exp(ln_pi_expected-0.5*trace(Bk*(R^-1))))/(exp(ln_pi_expected-0.5*trace(Bk*(R^-1)))+exp(ln_one_minus_pi_expected))
    for iz=1:dim
       
    z_t(iz)=1;
    pz1=exp( -.5*Bk(iz,iz)/R(iz,iz)+( psi(e_t(iz))-psi(e_t(iz)+f_t(iz))));
    z_t(iz)=0;
    pz0=exp( ( psi(f_t(iz))-psi(e_t(iz)+f_t(iz))) );
    rat=pz1/(pz1+pz0);
    
        if(rat>.5)
        z_t(iz)=1;
        else
        z_t(iz)=0;
        end
        
        if((z_t(iz))==0)
        Sig_zt(iz,iz)=inf;
        else
        Sig_zt(iz,iz)=R(iz,iz);
        end
    end   

    
    
    %% Updating e_t and f_t
    e_t = e_0 + z_t;
    f_t = f_0 + 1 - z_t; 
    
        
    d=norm(xp(:,j)-xp(:,j-1))/norm(xp(:,j-1));
    
    end
%     ind2(:,k)=z_t;
    z_t;
    xe=xp(:,j);
    Pe=Pp;    
    xpp=[xpp xe];
%     Ppp=[Ppp Pe];
end


    
    
%     rmse_ukf = sqrt(mean((x(1,:)-xpp(1,:)).^2));
% 
%     time_mse_sp=mean((x(1,:)-xpp(1,:)).^2+(x(3,:)-xpp(3,:)).^2);
    time_mse_sp=mean(mean((x(:,:)-xpp(:,:)).^2));
%     figure
%     plot(ind2)
    

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



function [inv_sig]=inver_sig_test(Sig)
    inv_sig=zeros(size(Sig));
    indc_inf=find(diag(Sig)~=inf);
    inv_sig(indc_inf,indc_inf)=Sig(indc_inf,indc_inf)^-1;
end




end