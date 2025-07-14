
% Implementation for the EMORF-II algorihtms
% Equations numbers refer to equations the paper
function xpp = EMORF_II(AA, B, a, bb, y, x, x_0, dt, Q, R)
    %—— Time-step and initial state setup ———————————————————————————
    N           = size(y, 2);    % Number of measurements
    xe          = x_0;           % Initial state estimate
    Pe          = Q;             % Initial state covariance
    

    %—— Initialize EMORF-II parameters ————————————————————————————————————
    theet = 0.5;                 
    alph  = a + 0.5;             

    %—— Unscented-transform weight setup ——————————————————————————
    n      = size(xe, 1);                            
    alpha  = 1;                                       
    beta   = 2;                                       
    kappa  = 0;                                       
    lambda = alpha^2 * (n + kappa) - n;              
    [WM, WC] = weights(n, alpha, beta, lambda);

    %—— Preallocate outputs and dimensions ————————————————————————
    xpp     = [];                
    dim     = size(y, 1);       


    %—— invert R once ——————————————————————              
    Rinv    = inv(R);           

                                  
for k = 1 : N
    % Initialize VB indicator vector
    I_t = ones(dim, 1);

    % Initialize beta vector
    for jjj = 1 : dim
        bet(jjj) = 1;    % initializing beta
    end

    % Implementing Eq.(28) for b_k
    E           = length((I_t ~= 1));
    sum_not_one = sum(I_t(I_t ~= 1));
    bb          = (E*a + AA - 1) / (sum_not_one + B);

    d = 1;

    %—— Prediction step ——————————————————————————————————
    % Generate sigma points
    SX = sigma_gen(xe, Pe, n, lambda);
    % Propagate through the dynamic model
    HX = dynamic_model(SX, dt);

    % Compute predicted mean and covariance (Eqs. 11 & 12)
    [xm, Pm] = mean_cov(WM, WC, HX);
    Pm = Pm + Q;

    j = 1;

    %—— Measurement prediction ————————————————————————————
    SX = sigma_gen(xm, Pm, n, lambda);
    HY = meas_model(SX, dim);

    %—— Update step ——————————————————————————————————————
    [mu, Uk] = mean_cov(WM, WC, HY);    
    C        = cross_cov(SX, HY, WC, xm, mu);

    Rkdi = Rinv;
    invS = Rkdi - Rkdi * ((eye(dim) + Uk*Rkdi)^-1) * Uk * Rkdi;
    K    = C * invS;

    % Implementing Eqs. (15) & (16)
    xp(:, j) = xm + K * (y(:, k) - mu);
    Pp       = Pm - C * K';

    % Compute mixture covariance W_k
    SX       = sigma_gen(xp(:, j), Pp, n, lambda);
    HY       = meas_model(SX, dim);
    [my, Yk] = mean_cov(WM, WC, HY);
    Bk       = (my - y(:, k)) * (my - y(:, k))' + Yk;  % this is W_k

    % Initialize measurement covariance for VB update
    Sig_It = R;

    %—— VB update per Eq.(23) —————————————————————————————
    for iz = 1 : dim
        F_4_1 = function4(Bk, I_t, R, iz);
        [R_iz_min, R_iz_iz_inv, det_R_iz_min_inv, trace_argument] = ...
            compute_matrix_operations(R, I_t, iz, Bk);

        alph    = a + 0.5;
        bet(iz) = bb + 0.5 * (Bk(iz, iz) / R(iz, iz));

        % Compute H and G (Eqs. 24 & 25)
        H = F_4_1 * theet;
        G = det_R_iz_min_inv * R_iz_iz_inv * ...
            exp(-0.5 * trace(trace_argument)) * ...
            (1 - theet) * ...
            ((bb^a / gamma(a)) / (bet(iz)^alph / gamma(alph)));

        % Update I_t per Eq.(23)
        if (H / G >= 1)
            I_t(iz) = 1;
        else
            I_t(iz) = (alph - 1) / bet(iz);
        end
    end

    %—— Update measurement covariance based on I_t —————————————
    for iz = 1 : dim
        if (I_t(iz) ~= 1)
            Sig_diag        = R(iz, iz) / ((alph - 1) / bet(iz));
            Sig_It(iz, :)   = 0;
            Sig_It(:, iz)   = 0;
            Sig_It(iz, iz)  = Sig_diag;
        end
    end

        % update b_k through Eq.(28)
        
        E = sum((I_t ~= 1));
        sum_not_one = sum(I_t(I_t ~= 1));
        bb = (E*a + AA - 1) / (sum_not_one + B);


    %---- VB loop --------------------------------------------------------
    while (d > 10^-4)
        j = j + 1;

        %---- update ------------------------------------------------------
        Rkdi = inver_sig_test_filcrem(Sig_It, I_t, eps);
        invS = Rkdi - Rkdi * ((eye(dim) + Uk * Rkdi)^-1) * Uk * Rkdi;
        K    = C * invS;
        % Eq.(15) & (16): state update inside the VB Loop
        xp(:, j) = xm + K * (y(:, k) - mu);
        Pp       = Pm - C * K';

        %---- Calculating W_k ----------------------------------------------
        SX       = sigma_gen(xp(:, j), Pp, n, lambda);
        HY       = meas_model(SX, dim);        
        [my, Yk] = mean_cov(WM, WC, HY);  
        Bk       = (my - y(:, k)) * (my - y(:, k))' + Yk;  % this is W_k

        %---- update I_t through Eq.(23) ----------------------------------
        for iz = 1 : dim
            F_4_1 = function4(Bk, I_t, R, iz);
            [R_iz_min, R_iz_iz_inv, det_R_iz_min_inv, trace_argument] = ...
                compute_matrix_operations(R, I_t, iz, Bk);

            alph    = a + 0.5;
            bet(iz) = bb + 0.5 * (Bk(iz, iz) / R(iz, iz));
            % Eqs.(24) & (25): compute H & G
            H = F_4_1 * theet;
            G = det_R_iz_min_inv * R_iz_iz_inv * ...
                exp(-0.5 * trace(trace_argument)) * ...
                (1 - theet) * ...
                (((bb^a) / gamma(a)) / ((bet(iz)^alph) / gamma(alph)));

            if (H / G >= 1)
                I_t(iz) = 1;
            else
                I_t(iz) = (alph - 1) / bet(iz);
            end
        end

        %---- update measurement covariance based on I_t ------------------
        Sig_It = R;  % reset base covariance
        for iz = 1 : dim
            if (I_t(iz) ~= 1)
                Sig_diag       = R(iz, iz) / ((alph - 1) / bet(iz));
                Sig_It(iz, :)  = 0;
                Sig_It(:, iz)  = 0;
                Sig_It(iz, iz) = Sig_diag;
            end
        end

        %---- update b_k through Eq.(28) ----------------------------------
        E           = sum((I_t ~= 1));
        sum_not_one = sum(I_t(I_t ~= 1));
        bb          = (E * a + AA - 1) / (sum_not_one + B);

        if (j > 20000)
            break
        end
   
    % checking convergence
        
    d=norm(xp(:,j)-xp(:,j-1))/norm(xp(:,j-1));
    
    end

    xe=xp(:,j);
    Pe=Pp;    
    xpp=[xpp xe];
%     Ppp=[Ppp Pe];
end


    time_mse_sp=mean(mean((x(:,:)-xpp(:,:)).^2));
    
end
%-----------helper functions--------------------
% this function is used for generating the sigma points
function SS=sigma_gen(x,P,n,lambda)  
    P=nearestSPD(P);
    A = chol(P,'lower');
    SS = [zeros(size(x)) A -A];
    SS = sqrt(n + lambda)*SS + repmat(x,1,size(SS,2));
end
% this function passes the sigma points through the state model
function HX=dynamic_model(SX,dt)
     HX=zeros(size(SX,1),size(SX,2));
        for i=1:size(SX,2)
        w_t=SX(5,i);
        F_t=[1 sin(w_t*dt)/w_t 0 (cos(w_t*dt)-1)/w_t 0;...
        0 cos(w_t*dt) 0 -sin(w_t*dt) 0;...
        0 (1-cos(w_t*dt))/w_t 1 sin(w_t*dt)/w_t 0;...
        0 sin(w_t*dt) 0 cos(w_t*dt) 0;...
        0 0 0 0 1];
        HX(:,i) = F_t*SX(:,i);
        end
end
% this function calculates the mean and variance for the state
function [x,P]=mean_cov(WM,WC,HX)
x=0;P=0;
        for i=1:size(HX,2)
            x = x + WM(i) * HX(:,i);
        end
        for i=1:size(HX,2)
            P = P + WC(i) * (HX(:,i) - x) * (HX(:,i) - x)';
        end
end
% calculates the cross covariance
function [C]=cross_cov(SX,HY,WC,xm,mu)
    C  = zeros(size(SX,1),size(HY,1));
    for i=1:size(SX,2)
            C = C + WC(i) * (SX(:,i) - xm) * (HY(:,i) - mu)';
    end
end
% passes sigma points through the measurement model
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
    
    for jt=1:2*n+1
        if jt==1
            wm = lambda / (n + lambda);
            wc = lambda / (n + lambda) + (1 - alpha^2 + beta);
        else
            wm = 1 / (2 * (n + lambda));
            wc = wm;
        end
        WM(jt) = wm;
        WC(jt) = wc;
    end
end
   

function [inv_sig] = inver_sig_test_filcrem(Sig,It,ep)
    inv_sig = zeros(size(Sig));
    
    % Find indices where It == 1
    indc_one = find(It == 1);
    
    % Find indices where It ~= 1
    indc_other = find(It ~= 1);
    
    % Full inverse for submatrix where It == 1
    if ~isempty(indc_one)
        inv_sig(indc_one, indc_one) = inv(Sig(indc_one, indc_one));
    end
    
    % Diagonal inverse for submatrix where It ~= 1
    if ~isempty(indc_other)
        inv_sig(indc_other, indc_other) = diag(1 ./ diag(Sig(indc_other, indc_other)));
    end
end

 function M = function4(Bk, I_ti, Ri, iz)
    % Get the dimension of R
    dimer = size(Ri, 1);

    % Initialize R_nominal as a copy of R
    R_nominal = Ri;

    % Modify R_nominal based on I_ti and iz
    for ii = 1:dimer
        if ii == iz
            % For the iz index, leave the iz-th row and column unchanged
            R_nominal(ii, ii) = Ri(ii, ii);  % Set specific index to Ri(ii, ii)
        elseif I_ti(ii) ~= 1
            % For indices where I_ti(i) != 1, set the i-th row and column to zero
            % except for the diagonal element
            R_nominal(ii, :) = 0;  % Set i-th row to zero
            R_nominal(:, ii) = 0;  % Set i-th column to zero
            R_nominal(ii, ii) = Ri(ii, ii) / I_ti(ii);  % Set diagonal element
        else
            % For indices where I_ti(i) == 1, leave the row and column unchanged
            R_nominal(ii, ii) = Ri(ii, ii)/ I_ti(ii);  % Keep diagonal element unchanged
        end
    end

    
    % Compute determinant and inverse of R_nominal
    det_R = det(R_nominal);
    %replacing the direct inversion with a somewhat more efficient
    %inversion
    %inv_Ri = inv(R_nominal);
    I_ti_copy = I_ti;
    I_ti_copy(iz) = 1;

    inv_Ri = inver_sig_test_filcrem(R_nominal,I_ti_copy,10); % this seems to be only optmization we can do here !
    
    % disp('New Way')
    % inv_Ri

    % Compute seg_A
    % this was just a single number there should not be a problem here
    seg_A = det_R^(-1/2);

    % Compute seg_B
    trace_term = trace(Bk * inv_Ri);
    seg_B = exp(-0.5 * trace_term);

    % Compute final matrix M
    M = seg_A * seg_B;
end


    function [R_iz_min, R_iz_iz_inv, det_R_iz_min_inv, trace_argument] = compute_matrix_operations(Ri, I_ti, izi, Bki)
    % Get the dimension of R
    dimer = size(Ri, 1);

    % Initialize R_nominal as a copy of R
    R_nominal = Ri;

    % Modify R_nominal based on I_ti and izi
    for ii = 1:dimer
        % if ii == izi
        %     % For the izi index, leave the izi-th row and column unchanged
        %     R_nominal(ii, ii) = Ri(ii, ii);  % Set diagonal element to Ri(ii, ii)
        %     R_nominal(ii,:) = 0;
        %     R_nominal(:,ii) = 0;
        if I_ti(ii) ~= 1
            % For indices where I_ti(i) != 1, set the i-th row and column to zero
            % except for the diagonal element
            R_nominal(ii, :) = 0;  % Set i-th row to zero
            R_nominal(:, ii) = 0;  % Set i-th column to zero
            R_nominal(ii, ii) = Ri(ii, ii) / I_ti(ii);  % Set diagonal element
        else
            % For indices where I_ti(i) == 1, leave the row and column unchanged
            R_nominal(ii, ii) = Ri(ii, ii)/ I_ti(ii);  % Keep diagonal element unchanged
        end
    end

    % Extract the iz-th diagonal entry from the original R
    R_iz_iz = Ri(izi, izi);
    
    % Compute R_iz_iz^(-1/2)
    % this is the inverting of a sngle number so no problem here !
    R_iz_iz_inv = R_iz_iz^(-1/2);

    % Remove the iz-th row and column from R_nominal to get R_iz_min
    R_iz_min = R_nominal;
    I_iz_min = I_ti;
    R_iz_min(izi, :) = [];  % Remove iz-th row
    R_iz_min(:, izi) = [];  % Remove iz-th column

    I_iz_min(izi) = []; %removing the iz-th column to go along with R_iz_min
    % disp('the izi value is')
    % izi

    % Compute det(R_iz_min)^(-1/2)
    % also the inversion of a single number 
    % which is the deteminent so no big deal

    det_R_iz_min = det(R_iz_min);
    det_R_iz_min_inv = det_R_iz_min^(-1/2);

    % Create a copy of Bk where the iz-th row and column are removed
    Bk_iz_min = Bki;
    Bk_iz_min(izi, :) = [];  % Remove iz-th row
    Bk_iz_min(:, izi) = [];  % Remove iz-th column

    % Compute trace_argument = Bk_iz_min * inv(R_iz_min)

    % we could technically imporive this inverson lets try

    % disp('New Style')
    % inver_sig_test_filcrem(R_iz_min,I_iz_min,10)


    R_iz_min_inv = inver_sig_test_filcrem(R_iz_min,I_iz_min,10);

    trace_argument = Bk_iz_min * R_iz_min_inv;
end

