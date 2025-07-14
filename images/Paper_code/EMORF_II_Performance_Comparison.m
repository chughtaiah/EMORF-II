%% Using new library!

clear all;
close all;

%% Parameters
T = 100;                        % Time steps
alpha = 1000;                   % Parameter for target tracking
dim = 5;                        % Dimension (fixed in this example)
num_runs = 10;                 % Number of MC runs for each lambda value
lambda_contam_values = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6];  % Contamination levels

%% Preallocate RMSE storage (num_runs x number of lambda values)
% Order of algorithms:
% 1. Ideal UKF
% 2. EMORF-II
% 3. EMORF
% 4. Gen. VBKF (N = 1)
% 5. Gen. VBKF (N = 10)
% 6. Ind. VBKF
nAlgo = 6;
rmse_ukf        = zeros(num_runs, length(lambda_contam_values));
rmse_emorf_asor = zeros(num_runs, length(lambda_contam_values));
rmse_emorf      = zeros(num_runs, length(lambda_contam_values));
rmse_vbkf_n1    = zeros(num_runs, length(lambda_contam_values));
rmse_vbkf_n10   = zeros(num_runs, length(lambda_contam_values));
rmse_vbkf_ind   = zeros(num_runs, length(lambda_contam_values));

%% Compute RMSE for each lambda value and each run
for lambda_idx = 1:length(lambda_contam_values)
    lambda_contam = lambda_contam_values(lambda_idx);
    for run = 1:num_runs
        % Generate data
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
        
        % Compute RMSE using all dimensions:
        rmse_ukf(run, lambda_idx)        = sqrt(mean((xpp_ukf - xout).^2, 'all'));
        rmse_emorf_asor(run, lambda_idx)   = sqrt(mean((xpp_emorf_a - xout).^2, 'all'));
        rmse_emorf(run, lambda_idx)        = sqrt(mean((xpp_emorf - xout).^2, 'all'));
        rmse_vbkf_n1(run, lambda_idx)      = sqrt(mean((vbkf_n_1 - xout).^2, 'all'));
        rmse_vbkf_n10(run, lambda_idx)     = sqrt(mean((vbkf_n_10 - xout).^2, 'all'));
        rmse_vbkf_ind(run, lambda_idx)     = sqrt(mean((vbkf_ind_xp - xout).^2, 'all'));
    end
end

%% Combine RMSE data into a 3D matrix
% Desired size: [num_runs x number of lambda values x number of algorithms]
nLambda = length(lambda_contam_values);
data3D = zeros(num_runs, nLambda, nAlgo);
for lambda_idx = 1:nLambda
    % Order: Ideal UKF, EMORF-II, EMORF, Gen. VBKF (N = 1), Gen. VBKF (N = 10), Ind. VBKF
    data3D(:, lambda_idx, :) = [ rmse_ukf(:, lambda_idx), ...
                                  rmse_emorf_asor(:, lambda_idx), ...
                                  rmse_emorf(:, lambda_idx), ...
                                  rmse_vbkf_n1(:, lambda_idx), ...
                                  rmse_vbkf_n10(:, lambda_idx), ...
                                  rmse_vbkf_ind(:, lambda_idx) ];
end

%% Font Control
% Adjust these variables to control various font sizes and axis limits.
titleFontSize   = 30;
legendFontSize  = 25;
xLabelFontSize  = 30;
yLabelFontSize  = 30;
xTickFontSize   = 30;
yTickFontSize   = 30;
ylimValue       = [0, 100];  % Set desired y-axis limits

%% Prepare Data for Box Plots
% Create a numeric array Y_box of size [num_runs*nAlgo x nLambda] where each column
% corresponds to one lambda (contamination level).
Y_box = [];
for j = 1:nLambda
    % For each lambda level, extract a [num_runs x nAlgo] slice,
    % then reshape it so that the data from each algorithm is stacked.
    Y_box = [Y_box, reshape(squeeze(data3D(:, j, :)), [], 1)];
end

% Create the grouping vector so that it matches the stacking order.
% This yields a vector where, for each lambda level, the first num_runs entries 
% are from algorithm 1, the next num_runs from algorithm 2, etc.
group_vec = reshape(repmat(1:nAlgo, num_runs, 1), [], 1);

%% Define Algorithm Names and Colors
algo_names = {'Ideal UKF', 'EMORF-II', 'EMORF', 'Gen. VBKF (N = 1)', 'Gen. VBKF (N = 10)', 'Ind. VBKF'};
% Define a color matrix with one row per algorithm.
% MATLAB's lines() returns a consistent set of colors each time.
c = lines(nAlgo);

% Define x-axis labels for the lambda contamination levels.
lambda_labels = arrayfun(@num2str, lambda_contam_values, 'UniformOutput', false);

figure('Name', 'RMSE Box Plots','WindowStyle','docked');
hBox = daboxplot(Y_box, 'groups', group_vec, ...
    'legend', algo_names, ...
    'xtlabels', lambda_labels, ...
    'fill', 0, ...         % Use non-filled boxplots (only outlines)
    'colors', c, ...       % Each row corresponds to one algorithm’s RGB color
    'whiskers', 1, ...     % Enable whiskers
    'scatter', 0);         % No scatter overlay

%% Turn on the box to draw a full rectangle around the plot
box on;

%% Add dashed vertical separator lines between conditions (without legends)
hold on;
nConditions = length(lambda_labels);
for j = 1:(nConditions-1)
    % Separator line position: between conditions j and j+1.
    x_sep = j + 0.5;
    % Draw a dashed line (LineStyle '--') across the entire y-axis.
    % 'HandleVisibility','off' prevents these lines from appearing in the legend.
    line([x_sep, x_sep], ylim, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

%% Apply Font and Axis Controls
xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', xLabelFontSize);
ylabel('RMSE', 'FontSize', yLabelFontSize);
set(gca, 'FontSize', xTickFontSize);  % Sets tick label font size
ylim(ylimValue);
if isfield(hBox, 'lg') && ~isempty(hBox.lg)
    legend('FontSize', legendFontSize, 'Location', 'northwest', 'NumColumns', 2);
end