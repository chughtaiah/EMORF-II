function xpp =robust_EMORF_self_modular_ind_sens_nsensors_xpp(y,x,x_0,dt,Q,R)


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
    dim=size(y,1);
    eps=10^-6;
    Rinv=inv(R);
    
for k=1:N
        d=1;
        %generate sigmas
        SX=sigma_gen(xe,Pe,n,lambda);        
        % Propagate through the dynamic model 
        HX=dynamic_model(SX,dt);        
        
        % Compute the predicted mean and covariance
        [xm,Pm]=mean_cov(WM,WC,HX);   
        Pm = Pm + Q;

        I_t=ones(dim,1);
        j=1;
    
        % Form sigma points for measurement step and
        % propagate throught the measurement model
        SX=sigma_gen(xm,Pm,n,lambda);        
        HY=meas_model(SX,dim);   
    
        
        [mu,Uk]=mean_cov(WM,WC,HY);    
        C=cross_cov(SX,HY,WC,xm,mu);
        
    Rkdi=Rinv;
    
    invS=Rkdi-Rkdi*((eye(dim)+Uk*Rkdi)^-1)*Uk*Rkdi;
    K=C*invS;

    
    xp(:,j)=xm+K*(y(:,k)-mu);
    Pp=Pm-C*K';

    SX=sigma_gen(xp(:,j),Pp,n,lambda);
    HY=meas_model(SX,dim);
        
    [my,Yk]=mean_cov(WM,WC,HY);              
    Bk=(my-y(:,k))*(my-y(:,k))'+Yk;
    
    Sig_It=R;
    
     for iz=1:dim
%     I_tm=I_t;
    
%     I_t1=I_t;    
%     I_t2=I_t;
%     
%     I_t1(iz)=1;
%     Sig_temp1=Sig_It;
%     Sig_temp1(iz,:)=R(iz,:);
%     Sig_temp1(:,iz)=R(:,iz);   
%     
% %     pz1=trace(Bk*inver_sig_test_filcrem(Sig_temp1,I_t1))+log(det(Sig_temp1));
% 
% 
%     I_t2(iz)=0;
%     Sig_temp2=Sig_It;
%     Sig_temp2(iz,:)=0;
%     Sig_temp2(:,iz)=0;
%     Sig_temp2(iz,iz)=R(iz,iz)/eps;


        ai=1:dim;
        It_temp=I_t;
        It_temp(iz)=[];
        ai(iz)=[];
%         Sigg=Sig_temp1(ai,ai);
%         b=Sig_temp1(iz,ai);
%         dif=trace(Bk*(inver_sig_test_filcrem(Sig_temp1,I_t1)-inver_sig_test_filcrem(Sig_temp2,I_t2)))+log(eps)+log(det(Sigg-b'*b/(R(iz,iz))))-log(det(Sigg))
 
        Aa=(R(iz,iz));
        Bb=R(iz,ai);
        Cc=Bb';
        Dd=Sig_It(ai,ai); % isnt this incorrect,  shouldnt this be without the  iz dimension ?
%         [Aa Bb;Cc Dd];
        in_Dd=inver_sig_test_filcrem(Dd,It_temp);
        ABDC=1/(Aa-Bb*in_Dd*Cc);
        Mat11=ABDC-eps/Aa; 
        Mat12=-ABDC*Bb*in_Dd;
        Mat21=Mat12';
        Mat22=in_Dd*Cc*ABDC*Bb*in_Dd;
        Mat_unswap=zeros(dim,dim);
%         Mat_swap=[Mat11 Mat12; Mat21 Mat22];
        Mat_unswap(iz,iz)=Mat11;
        Mat_unswap(ai,ai)=Mat22;
        Mat_unswap(iz,ai)=Mat12;
        Mat_unswap(ai,iz)=Mat21;
%         trace(Bk*Mat_unswap);
%         trace(Bk*(inver_sig_test_filcrem(Sig_temp1,I_t1)-inver_sig_test_filcrem(Sig_temp2,I_t2)))
%         dif=trace(Bk*(inver_sig_test_filcrem(Sig_temp1,I_t1)-inver_sig_test_filcrem(Sig_temp2,I_t2)))+log(eps)+log(det(eye(dim-1)-b'*b*inver_sig_test_filcrem(Sigg,It_temp)/(R(iz,iz))));
        dif=trace(Bk*(Mat_unswap))+log(eps)+log(det(eye(dim-1)-Bb'*Bb*in_Dd/(R(iz,iz))));


        if(dif<0)   
        I_t(iz)=1;
        else
        I_t(iz)=0;
        end
        
        if((I_t(iz))==0)
        Sig_diag=R(iz,iz)/eps;
        Sig_It(iz,:)=0;
        Sig_It(:,iz)=0;
        Sig_It(iz,iz)=Sig_diag;
        else
        Sig_It(iz,:)=R(iz,:);
        Sig_It(:,iz)=R(:,iz);
        end
        
        
    end 

%         
%     end
    
    %% Updating e_t and f_t


    
    while(d>10^-4)
    j=j+1;

    Rkdi=inver_sig_test_filcrem(Sig_It,I_t);
    invS=Rkdi-Rkdi*((eye(dim)+Uk*Rkdi)^-1)*Uk*Rkdi;
    K=C*invS;
    
    xp(:,j)=xm+K*(y(:,k)-mu);
    Pp=Pm-C*K';    

    SX=sigma_gen(xp(:,j),Pp,n,lambda);
    HY=meas_model(SX,dim);        
    [my,Yk]=mean_cov(WM,WC,HY);  
        
    
    Bk=(my-y(:,k))*(my-y(:,k))'+Yk;

%     d_i=1;
%     
%     while d_i>10^-4
    
    for iz=1:dim


        ai=1:dim;
        It_temp=I_t;
        It_temp(iz)=[];
        ai(iz)=[];

        Aa=(R(iz,iz));
        Bb=R(iz,ai);
        Cc=Bb';
        Dd=Sig_It(ai,ai);

        in_Dd=inver_sig_test_filcrem(Dd,It_temp);
        ABDC=1/(Aa-Bb*in_Dd*Cc);
        Mat11=ABDC-eps/Aa; 
        Mat12=-ABDC*Bb*in_Dd;
        Mat21=Mat12';
        Mat22=in_Dd*Cc*ABDC*Bb*in_Dd;
        Mat_unswap=zeros(dim,dim);

        Mat_unswap(iz,iz)=Mat11;
        Mat_unswap(ai,ai)=Mat22;
        Mat_unswap(iz,ai)=Mat12;
        Mat_unswap(ai,iz)=Mat21;

        dif=trace(Bk*(Mat_unswap))+log(eps)+log(det(eye(dim-1)-Bb'*Bb*in_Dd/(R(iz,iz))));


        if(dif<0)   
        I_t(iz)=1;
        else
        I_t(iz)=0;
        end
        
        if((I_t(iz))==0)
        Sig_diag=R(iz,iz)/eps;
        Sig_It(iz,:)=0;
        Sig_It(:,iz)=0;
        Sig_It(iz,iz)=Sig_diag;
        else
        Sig_It(iz,:)=R(iz,:);
        Sig_It(:,iz)=R(:,iz);
        end
        
        
    end 
%         d_i=norm(I_t-I_tm)/norm(I_tm);
        
%     end
    
    if(j>20000)
        break
    end
    
    %% Updating e_t and f_t
        
    d=norm(xp(:,j)-xp(:,j-1))/norm(xp(:,j-1));
    
    end
%     ind2(:,k)=z_t;
    I_t;
    xe=xp(:,j);
    Pe=Pp;    
    xpp=[xpp xe];
%     Ppp=[Ppp Pe];
end


%     figure
%     subplot(2,1,1)
%     plot(x(1,:),'r-');
%     hold on 
%     plot(xpp(1,:));
%     title('UKF estimate');
%     
%     
%     subplot(2,1,2)
%     plot(x(3,:),'r-');
%     hold on 
%     plot(xpp(3,:));
    
%     rmse_ukf = sqrt(mean((x(1,:)-xpp(1,:)).^2));
% %     size(x(1,:))
% %     size(xpp(1,:))
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
    

function [inv_sig]=inver_sig_test_filcrem(Sig,It)
    inv_sig=zeros(size(Sig));
    indc_one=find(It==1);
    indc_z=find(It==0);
%     indc=mod(indc,size(Sig,1));
    if ~isempty(indc_one)
    inv_sig(indc_one,indc_one)=inv(Sig(indc_one,indc_one));
    end
    
    if ~isempty(indc_z)
    inv_sig(indc_z,indc_z)=diag(1./diag(Sig(indc_z,indc_z)));
    end
%      inv_sig(indc_z,indc_z)=0;
end




end