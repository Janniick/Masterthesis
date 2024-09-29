function mutual_information_band = mutual_information_band(delta_ts, theta_ts, alpha_ts, beta_ts, spindles_ts)
    % Function to calculate mutual information between different frequency bands.
    %
    % Inputs:
    %   - delta_ts: time series for delta band
    %   - theta_ts: time series for theta band
    %   - alpha_ts: time series for alpha band
    %   - beta_ts: time series for beta band
    %   - spindles_ts: time series for spindles band
    %
    % Output:
    %   - mutual_information_band: a struct containing mutual information between the various frequency bands.

    % Initialize struct to store mutual information
    mutual_information_band = struct();

    % Calculate mutual information between all band combinations
    mutual_information_band.delta_theta = gcmi_cc(delta_ts(:), theta_ts(:));
    mutual_information_band.delta_alpha = gcmi_cc(delta_ts(:), alpha_ts(:));
    mutual_information_band.delta_beta = gcmi_cc(delta_ts(:), beta_ts(:));
    mutual_information_band.delta_spindles = gcmi_cc(delta_ts(:), spindles_ts(:));
    
    mutual_information_band.theta_alpha = gcmi_cc(theta_ts(:), alpha_ts(:));
    mutual_information_band.theta_beta = gcmi_cc(theta_ts(:), beta_ts(:));
    mutual_information_band.theta_spindles = gcmi_cc(theta_ts(:), spindles_ts(:));
    
    mutual_information_band.alpha_beta = gcmi_cc(alpha_ts(:), beta_ts(:));
    mutual_information_band.alpha_spindles = gcmi_cc(alpha_ts(:), spindles_ts(:));
    
    mutual_information_band.beta_spindles = gcmi_cc(beta_ts(:), spindles_ts(:));

end
