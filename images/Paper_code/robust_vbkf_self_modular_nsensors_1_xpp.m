function xpp=robust_vbkf_self_modular_nsensors_1_xpp(y,x,x_0,dt,Q,R)


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
    % ratio=zeros(dim,1);

    MC_N=1;   
    Sig_zt=zeros(dim,dim,MC_N);

for k=1:N
        d=1;
        %generate sigmas
        SX=sigma_gen(xe,Pe,n,lambda);        
        % Propagate through the dynamic model 
        HX=dynamic_model(SX,dt);        
        
        % Compute the predicted mean and covariance
        [xm,Pm]=mean_cov(WM,WC,HX);   
        Pm = Pm + Q;

        z_t=ones(dim,MC_N);
        e_t = e_0*ones(dim,1);
        e_t = e_t +ones(dim,1);
        f_t = f_0*ones(dim,1);
    
        j=1;
    
        % Form sigma points for measurement step and
        % propagate throught the measurement model
        SX=sigma_gen(xm,Pm,n,lambda);        
        HY=meas_model(SX,dim);   
    
        
        [mu,Uk]=mean_cov(WM,WC,HY);    
        C=cross_cov(SX,HY,WC,xm,mu);
        
    
    for mc=1:MC_N
        Sig_zt(:,:,mc)=R;
    end
    K = C*((Uk + R)^-1);  
    
    % K = C*((Uk + Sig_zt)^-1);  

    
    xp(:,j)=xm+K*(y(:,k)-mu);
    Pp=Pm-C*K';

    SX=sigma_gen(xp(:,j),Pp,n,lambda);
    HY=meas_model(SX,dim);
        
    [my,Yk]=mean_cov(WM,WC,HY);              
    Bk=(my-y(:,k))*(my-y(:,k))'+Yk;
    
    
    for iz=1:dim

    [Bk_mi_mi,Bk_i_mi,Bk_mi_i,Bk_i_i]=Sub(Bk,iz);
    
    E_tri_t1=0;
    E_tri_t2=0;
    E_tri_t3=0;
    E_lam_t4=0;

    for mc=1:MC_N
        [t1,t2,t3,t4]=E_tri(Sig_zt(:,:,mc),iz);
        E_tri_t1=E_tri_t1+t1;
        E_tri_t2=E_tri_t2+t2;
        E_tri_t3=E_tri_t3+t3;
        E_lam_t4=E_lam_t4+t4;
    end

    E_tri_t1=E_tri_t1/MC_N;
    E_tri_t2=E_tri_t2/MC_N;
    E_tri_t3=E_tri_t3/MC_N;
    E_lam_t4=E_lam_t4/MC_N;

    z_t(iz,:)=1;
    log_deter=0;
    
    for mc=1:MC_N
        % Sig_zt(iz,iz,mc)=R(iz,iz);
        % Sig_zt(iz,:,mc)=R(iz,:);
        % Sig_zt(:,iz,mc)=R(:,iz);
        log_deter=log_deter+log(det_C(Sig_zt(:,:,mc),z_t(:,mc)));
    end
    log_deter=log_deter/MC_N;

    % pz11=-.5*trace(Bk_mi_mi*E_tri_t1)+.5*trace(Bk_i_mi*E_tri_t2)+.5*trace(Bk_i_mi*E_tri_t3)-.5*trace(Bk_i_i*E_lam_t4)+( psi(e_t(iz))-psi(e_t(iz)+f_t(iz)) )-.5*log(2*pi)-.5*( log_deter );
    pz1=exp( -.5*trace(Bk_mi_mi*E_tri_t1)+.5*trace(Bk_i_mi*E_tri_t2)+.5*trace(Bk_i_mi*E_tri_t3)-.5*trace(Bk_i_i*E_lam_t4)+( psi(e_t(iz))-psi(e_t(iz)+f_t(iz)) )-.5*log(2*pi)-.5*( log_deter ) );
    z_t(iz,:)=0;
    log_deter=0;
    for mc=1:MC_N
        % Sig_zt(iz,iz,mc)=inf;
        % Sig_zt(iz,:,mc)=0;
        % Sig_zt(:,iz,mc)=0;
        % Sig_zt(iz,iz,mc)=inf;

        log_deter=log_deter+log(det_C(Sig_zt(:,:,mc),z_t(:,mc)));
    end
    log_deter=log_deter/MC_N;
    % pz00=( psi(f_t(iz))-psi(e_t(iz)+f_t(iz)) )-.5*( log_deter ) ;
    pz0=exp( ( psi(f_t(iz))-psi(e_t(iz)+f_t(iz)) )-.5*( log_deter ) );
    % rat_exp=1/(1+exp(pz00-pz11))
    rat=pz1/(pz1+pz0);
    
    % if(rat>.5)
    %     rat=1;
    % end
    % ratio(iz,1)=rat;



    % if(rat>.5)
    %     z_t(iz)=1;
    %     % Sig_zt(iz,iz)=R(iz,iz);
    %     Sig_zt(iz,:)=R(iz,:);
    %     Sig_zt(:,iz)=R(:,iz);
    %     % Sig_It(iz,:)=0;
    %     % Sig_It(:,iz)=0;
    %     else
    %     z_t(iz)=0;
    %     Sig_zt(iz,:)=0;
    %     Sig_zt(:,iz)=0;
    %     Sig_zt(iz,iz)=inf;
    % end

    for mc=1:MC_N
        if(rat>rand)
            z_t(iz,mc)=1;
            % Sig_zt(iz,:,mc)=R(iz,:);
            % Sig_zt(:,iz,mc)=R(:,iz);
        else
            z_t(iz,mc)=0;
            Sig_zt(iz,:,mc)=0;
            Sig_zt(:,iz,mc)=0;
            Sig_zt(iz,iz,mc)=inf;
        end

    end
    e_t(iz,1) = e_0 + rat;
    f_t(iz,1) = f_0 + 1 - rat;
 
    end 

    
    
    %% Updating e_t and f_t
    % e_t = e_0 + z_t;
    % f_t = f_0 + 1 - z_t; 
    % ratio=mean(z_t,2);
    % 
    % e_t = e_0 + ratio;
    % f_t = f_0 + 1 - ratio; 

    
    while(d>10^-4)
    j=j+1;
    % Rkdi=inver_sig_test(Sig_zt);
    Rkdi=0;
    for mc=1:MC_N
      Rkdi=Rkdi+inver_sig_test(Sig_zt(:,:,mc));
    end
    Rkdi=Rkdi/MC_N;
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
   
    [Bk_mi_mi,Bk_i_mi,Bk_mi_i,Bk_i_i]=Sub(Bk,iz);
    
    E_tri_t1=0;
    E_tri_t2=0;
    E_tri_t3=0;
    E_lam_t4=0;

    for mc=1:MC_N
    [t1,t2,t3,t4]=E_tri(Sig_zt(:,:,mc),iz);
    E_tri_t1=E_tri_t1+t1;
    E_tri_t2=E_tri_t2+t2;
    E_tri_t3=E_tri_t3+t3;
    E_lam_t4=E_lam_t4+t4;
    end

    E_tri_t1=E_tri_t1/MC_N;
    E_tri_t2=E_tri_t2/MC_N;
    E_tri_t3=E_tri_t3/MC_N;
    E_lam_t4=E_lam_t4/MC_N;

    z_t(iz,:)=1;
    

    log_deter=0;
    
    for mc=1:MC_N
        % Sig_zt(iz,iz,mc)=R(iz,iz);
        % Sig_zt(iz,:,mc)=R(iz,:);
        % Sig_zt(:,iz,mc)=R(:,iz);

        log_deter=log_deter+log(det_C(Sig_zt(:,:,mc),z_t(:,mc)));
    end
    log_deter=log_deter/MC_N;

    % pz11=-.5*trace(Bk_mi_mi*E_tri_t1)+.5*trace(Bk_i_mi*E_tri_t2)+.5*trace(Bk_i_mi*E_tri_t3)-.5*trace(Bk_i_i*E_lam_t4)+( psi(e_t(iz))-psi(e_t(iz)+f_t(iz)) )-.5*log(2*pi)-.5*( log_deter );
    pz1=exp( -.5*trace(Bk_mi_mi*E_tri_t1)+.5*trace(Bk_i_mi*E_tri_t2)+.5*trace(Bk_i_mi*E_tri_t3)-.5*trace(Bk_i_i*E_lam_t4)+( psi(e_t(iz))-psi(e_t(iz)+f_t(iz)) )-.5*log(2*pi)-.5*( log_deter ) );
    
    z_t(iz,:)=0;

    

    log_deter=0;
    for mc=1:MC_N
        % Sig_zt(iz,:,mc)=0;
        % Sig_zt(:,iz,mc)=0;
        % Sig_zt(iz,iz,mc)=inf;
        log_deter=log_deter+log(det_C(Sig_zt(:,:,mc),z_t(:,mc)));
    end
    log_deter=log_deter/MC_N;
    % pz00=( psi(f_t(iz))-psi(e_t(iz)+f_t(iz)) )-.5*( log_deter ) ;
    pz0=exp( ( psi(f_t(iz))-psi(e_t(iz)+f_t(iz)) )-.5*( log_deter ) );
    rat=pz1/(pz1+pz0);
    % pz00
    % pz11
    % pz00-pz11
    % exp(pz00-pz11)
    % rat_exp=1/(1+exp(pz00-pz11))
    

    % 
    % for mc=1:MC_N
    %     if(rat>rand)
    %         z_t(iz,mc)=1;
    %         Sig_zt(iz,iz,mc)=R(iz,iz);
    %     else
    %         z_t(iz,mc)=0;
    %         Sig_zt(iz,iz,mc)=inf;
    %     end
    % 
    % end
    
    for mc=1:MC_N
        if(rat>rand)
            z_t(iz,mc)=1;
            % Sig_zt(iz,:,mc)=R(iz,:);
            % Sig_zt(:,iz,mc)=R(:,iz);
            else
            z_t(iz,mc)=0;
            Sig_zt(iz,:,mc)=0;
            Sig_zt(:,iz,mc)=0;
            Sig_zt(iz,iz,mc)=inf;
        end
    end


    % ratio=mean(z_t,2);
    e_t(iz,1) = e_0 + rat;
    f_t(iz,1) = f_0 + 1 - rat;
    
    end  

    
    
    %% Updating e_t and f_t
    % e_t = e_0 + z_t;
    % f_t = f_0 + 1 - z_t; 
    % ratio=mean(z_t,2);
    % e_t = e_0 + ratio;
    % f_t = f_0 + 1 - ratio;
    % 
        
    d=norm(xp(:,j)-xp(:,j-1))/norm(xp(:,j-1));
    
    end
    % ratio;
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
    
function [A_mi_mi_s,A_i_mi_s,A_mi_i_s,A_i_i_s]=Sub(A,iii)

A_mi_mi_s=A;
A_mi_mi_s(iii,:)=[];
A_mi_mi_s(:,iii)=[];

A_i_mi_s=A(iii,:);
A_i_mi_s(iii)=[];

A_mi_i_s=A(:,iii);
A_mi_i_s(iii)=[];

A_i_i_s=A(iii,iii);
end

function [E_tri_t1,E_tri_t2,E_tri_t3,E_lam_t4]=E_tri(Sig,ii)
    
[Sig_mi_mi_s,Sig_i_mi_s,Sig_mi_i_s,Sig_i_i_s]=Sub(Sig,ii);

inver_mi_mi=inver_sig_test(Sig_mi_mi_s);
E_tri_t1=( inver_mi_mi*(Sig_mi_i_s*Sig_i_mi_s)*inver_mi_mi )/(Sig_i_i_s-Sig_i_mi_s*inver_mi_mi*Sig_mi_i_s);
E_tri_t2=2/(Sig_i_i_s)*inver_mi_mi*(Sig_mi_i_s) ;
E_tri_t3=( 2*inver_mi_mi*(Sig_mi_i_s*Sig_i_mi_s)*inver_mi_mi*Sig_mi_i_s )/(Sig_i_i_s^2-Sig_i_i_s*Sig_i_mi_s*inver_mi_mi*Sig_mi_i_s);
E_lam_t4=1/((Sig_i_i_s-Sig_i_mi_s*inver_mi_mi*Sig_mi_i_s));

end

function [inv_sig]=inver_sig_test(Sig)
    inv_sig=zeros(size(Sig));
    indc_inf=find(diag(Sig)~=inf);
    inv_sig(indc_inf,indc_inf)=Sig(indc_inf,indc_inf)^-1;
end

function detC=det_C(Sig,z_in)
    index=find(z_in==0);  
    Sig(index,:)=[];
    Sig(:,index)=[];  
    detC=det(Sig); 
end


end