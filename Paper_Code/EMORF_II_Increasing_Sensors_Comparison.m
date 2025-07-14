%% the better box plot code 
clear all;
close all;

%% Parameters
T = 100;                        % Time steps
alpha = 1000;                   % Parameter for target tracking
num_runs = 100;                  % Number of runs for each case
lambda_contam = 0.4;            % Fixed error probability
dims = [5,10,15,20];  % Sensor dimensions to test

%% Preallocate RMSE storage for each sensor dimension and algorithm
% Order of algorithms:
% 1. Ideal UKF
% 2. EMORF-II
% 3. EMORF
% 4. Gen. VBKF (N = 1)
% 5. Gen. VBKF (N = 10)
% 6. Ind. VBKF
nDims = length(dims);
nAlgo = 6;
rmse_ukf_all        = zeros(num_runs, nDims);
rmse_emorf_asor_all = zeros(num_runs, nDims);
rmse_emorf_all      = zeros(num_runs, nDims);
rmse_vbkf_n1_all    = zeros(num_runs, nDims);
rmse_vbkf_n10_all   = zeros(num_runs, nDims);
rmse_vbkf_ind_all   = zeros(num_runs, nDims);

%% Simulation: Loop over sensor dimensions
for d_idx = 1:nDims
    dim = dims(d_idx);
    
    % Preallocate for current dimension
    rmse_ukf        = zeros(num_runs, 1);
    rmse_emorf_asor = zeros(num_runs, 1);
    rmse_emorf      = zeros(num_runs, 1);
    rmse_vbkf_n1    = zeros(num_runs, 1);
    rmse_vbkf_n10   = zeros(num_runs, 1);
    rmse_vbkf_ind   = zeros(num_runs, 1);
    
    for run = 1:num_runs
        % Generate data for current dimension and fixed lambda_contam
        [R_out, Q_out, I, sig] = gen(dim, lambda_contam, T);
        [y_out, xout, xout_0, dt_out] = target_tracking_var_dim_nsensors(dim, R_out, Q_out, I, T, alpha, sig);
        
        % Run algorithms
        % EMORF
        xpp_emorf    = robust_EMORF_self_modular_ind_sens_nsensors_xpp(y_out, xout, xout_0, dt_out, Q_out, R_out);
        % EMORF-II
        AA = 10000; B = 1000; a = 1; b = 5000; theet = 0.5;
        xpp_emorf_a  = EMORF_II(AA, B, a, b, y_out, xout, xout_0, dt_out, Q_out, R_out);
        % Ideal UKF
        xpp_ukf      = ukf_ideal_self_modular_nsensors_xpp(y_out, xout, xout_0, dt_out, Q_out, R_out, I);
        % Gen. VBKF (N = 1)
        vbkf_n_1     = robust_vbkf_self_modular_nsensors_1_xpp(y_out, xout, xout_0, dt_out, Q_out, R_out);
        % Gen. VBKF (N = 10)
        vbkf_n_10    = robust_vbkf_self_modular_nsensors_10_xpp(y_out, xout, xout_0, dt_out, Q_out, R_out);
        % Ind. VBKF
        vbkf_ind_xp  = robust_vbkf_ind_self_modular_ind_sens_nsensors_xpp(y_out, xout, xout_0, dt_out, Q_out, R_out);
        
        % Compute RMSE (using the 1st and 3rd dimensions as an example)
        rmse_emorf(run)      = sqrt(mean((xpp_emorf(1,:) - xout(1,:)).^2 + (xpp_emorf(3,:) - xout(3,:)).^2, 'all'));
        rmse_emorf_asor(run) = sqrt(mean((xpp_emorf_a(1,:) - xout(1,:)).^2 + (xpp_emorf_a(3,:) - xout(3,:)).^2, 'all'));
        rmse_ukf(run)        = sqrt(mean((xpp_ukf(1,:) - xout(1,:)).^2 + (xpp_ukf(3,:) - xout(3,:)).^2, 'all'));
        rmse_vbkf_n1(run)    = sqrt(mean((vbkf_n_1(1,:) - xout(1,:)).^2 + (vbkf_n_1(3,:) - xout(3,:)).^2, 'all'));
        rmse_vbkf_n10(run)   = sqrt(mean((vbkf_n_10(1,:) - xout(1,:)).^2 + (vbkf_n_10(3,:) - xout(3,:)).^2, 'all'));
        rmse_vbkf_ind(run)   = sqrt(mean((vbkf_ind_xp(1,:) - xout(1,:)).^2 + (vbkf_ind_xp(3,:) - xout(3,:)).^2, 'all'));
    end
    
    % Store the RMSE for current dimension
    rmse_ukf_all(:, d_idx)        = rmse_ukf;
    rmse_emorf_asor_all(:, d_idx) = rmse_emorf_asor;
    rmse_emorf_all(:, d_idx)      = rmse_emorf;
    rmse_vbkf_n1_all(:, d_idx)    = rmse_vbkf_n1;
    rmse_vbkf_n10_all(:, d_idx)   = rmse_vbkf_n10;
    rmse_vbkf_ind_all(:, d_idx)   = rmse_vbkf_ind;
end

%% Combine RMSE data into a 3D matrix
% Desired size: [num_runs x number of sensor dimensions x number of algorithms]
data3D = zeros(num_runs, nDims, nAlgo);
for d_idx = 1:nDims
    data3D(:, d_idx, :) = [ rmse_ukf_all(:, d_idx), ...
                            rmse_emorf_asor_all(:, d_idx), ...
                            rmse_emorf_all(:, d_idx), ...
                            rmse_vbkf_n1_all(:, d_idx), ...
                            rmse_vbkf_n10_all(:, d_idx), ...
                            rmse_vbkf_ind_all(:, d_idx) ];
end

%% Font Control
titleFontSize   = 30;
legendFontSize  = 25;
xLabelFontSize  = 30;
yLabelFontSize  = 30;
xTickFontSize   = 30;
yTickFontSize   = 30;
ylimValue       = [0, 130];

%% Prepare Data for Box Plots
% Create a numeric array Y_box of size [num_runs*nAlgo x nDims] where each column
% corresponds to one sensor dimension.
Y_box = [];
for d_idx = 1:nDims
    Y_box = [Y_box, reshape(squeeze(data3D(:, d_idx, :)), [], 1)];
end

% Create the grouping vector so that it matches the stacking order.
% This yields a vector where, for each sensor dimension, the first num_runs entries 
% are from algorithm 1, the next num_runs from algorithm 2, etc.
group_vec = reshape(repmat(1:nAlgo, num_runs, 1), [], 1);

%% Define Labels, Algorithm Names, and Colors
% X-axis labels for sensor dimensions (converted to strings)
dim_labels = arrayfun(@num2str, dims, 'UniformOutput', false);

% Algorithm names (groups)
algo_names = {'Ideal UKF', 'EMORF-II', 'EMORF', 'Gen. VBKF (N = 1)', 'Gen. VBKF (N = 10)', 'Ind. VBKF'};

% Define a color matrix with one row per algorithm.
c = lines(nAlgo);

%% Create the Box Plot Figure Using daboxplot
figure('Name', 'RMSE Box Plots vs. Sensor Dimensions','WindowStyle','docked');
hBox = daboxplot(Y_box, 'groups', group_vec, ...
    'legend', algo_names, ...        % Legend based on algorithm names
    'xtlabels', dim_labels, ...        % X-axis labels for sensor dimensions
    'fill', 0, ...                   % Use non-filled boxplots (outlines only)
    'colors', c, ...                 % One row per algorithm
    'whiskers', 1, ...               % Enable whiskers
    'scatter', 0);                   % No scatter overlay

% Turn on the box to create a full border around the plot
box on;

%% Add dashed vertical separator lines between sensor dimensions (without legends)
hold on;
nConditions = length(dim_labels);
for j = 1:(nConditions-1)
    x_sep = j + 0.5;
    line([x_sep, x_sep], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

%% Apply Font and Axis Controls
xlabel('$m$', 'FontSize', xLabelFontSize, 'Interpreter', 'latex');
ylabel('RMSE', 'FontSize', yLabelFontSize);
set(gca, 'FontSize', xTickFontSize);
ylim(ylimValue);
if isfield(hBox, 'lg') && ~isempty(hBox.lg)
    legend('FontSize', legendFontSize, 'Location', 'northwest', 'NumColumns', 2);
end
