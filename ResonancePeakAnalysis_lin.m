function [Power_fit_store,lambda_fit_store,fitting,extract,param,resnorm] = ResonancePeakAnalysis_lin(data_power,data_lambda,FSR)
%%% This function fits all the resonance peaks in a measured MRR spectrum.
%%%
%%% For the fitting, a small portion of the spectrum is cut out around the
%%% resonance peak. A lorentz oscillator model is then used to fit the peak
%%% within this spectrum section where the spectrum is in units of W, i.e.
%%% the linearl spectrum is fitted such that the model parameters all have
%%% physical interpretations. The process is repeated for all peaks.
%%% This function fits the peaks in order from high to low wavelengths.
%%%
%%% Inputs:
%%% data_power = spectrum for specific device; array [mW]
%%% data_lambda = corresponding wavelength range; array [nm]
%%% FSR = approximated free spectral range; scalar [nm]
%%% 
%%% Outputs:
%%% Power_fit_store = measured power in fit range; matrix [mW]
%%% lambda_fit_store = wavelength in fit range; matrix [nm]
%%% fitting = fitted power in fit range; matrix [mW]
%%%     above data are stored for each peak; array row per peak
%%% extract = resonator data; cell containing
%%%     lambda_res = extracted resonance wavelength from fit [nm]
%%%     FWHM = extracted full-width half-max from fit [nm]
%%%     extinction = extinction ratio from measurement data [dB]
%%%     FSR = lambda_res spacings from fit [nm]
%%% param = fit parameters for the lorentz oscillator model
%%%
%%% Created by Matthias Vermeer, 5.23

%%% Adjusting x and y axis of spectrum, data{1} is measurement 1
PowerdB = real(double(10.*log10(data_power./1e-3)));
lambda = data_lambda;
dlambda = (max(lambda)-min(lambda))/(length(lambda)-1);

%%% Power smoothening
BW_smooth = 0.5;                    % BW of smooth filter [nm]
span = BW_smooth/dlambda;           % span for smooth filter
PowerSmooth = smooth(PowerdB,span);

%%% Defining leveling and fitting bandwidth
BWindex_level = round(FSR/dlambda); % interval for leveling of region 1
BW = FSR/4;                         % resonance fit BW
BWside = round(BW/2/dlambda);

%%% Leveling of data in FSR range
PowerSmooth = smooth(PowerdB,span);
slope = (PowerSmooth(1)-PowerSmooth(BWindex_level))/(lambda(1)-lambda(BWindex_level));
Power_level_loc = PowerdB(1:BWindex_level)-slope*(lambda(1:BWindex_level)-lambda(BWindex_level))-max(PowerdB(1:BWindex_level));
lambda_level_loc = lambda(1:BWindex_level);
Power_level = PowerdB-slope*(lambda-lambda(BWindex_level))-max(PowerdB(1:BWindex_level));
    
%%% Resonance peak fitting
%%% Determine "approx" resonance wavelength; as minimum
lambda_res_guess = sum(lambda_level_loc.*((Power_level_loc-min(Power_level_loc)).^2<1e-12));
lambda_res_index = sum((1:length(lambda)).*((lambda-lambda_res_guess).^2<1e-12));
%%% fit around "approx" resonance with Lorentzian distribution
if lambda_res_index>BWside  %%% check if the peak is not "too left"
Power_fit_dB = Power_level(lambda_res_index-BWside:lambda_res_index+BWside);
Power_fit = 10.^(Power_fit_dB/10)*1;
lambda_fit = lambda(lambda_res_index-BWside:lambda_res_index+BWside);
else    %%% otherwise, use same fitting bandwidth from index 1 on.
Power_fit_dB = Power_level(1:BWside*2+1);
Power_fit = 10.^(Power_fit_dB/10)*1;
lambda_fit = lambda(1:BWside*2+1); 
end
%%% lorentzian: Y(X) = P1./((X - P2).^2 + P3) + C.
P01 = -0.001;                % init: -1
P02 = lambda_res_guess;     
P03 = 0.001;                % init: 0.04
C0 = 1;                     % init: 40
P_initial = [P01 P02 P03 C0];
peak = 1;
[fitting(:,peak), param(:,peak), resnorm(:,peak)] = lorentzfit(lambda_fit,Power_fit,P_initial);
extinction_max(peak) = max(smooth(Power_fit_dB,span));
extinction_min(peak) = min(Power_fit_dB);
extinction_ER(peak) = max(smooth(Power_fit_dB,span))-min(Power_fit_dB);

Power_fit_store(peak,:) = Power_fit;
lambda_fit_store(peak,:) = lambda_fit;

while lambda_res_guess<(max(lambda)-FSR)
    index_peak = lambda_res_index + BWindex_level;
    peak = peak + 1;

    %%% Leveling of data in FSR range
    if (length(lambda)-index_peak)>BWside  %%% check if the peak is not "too right"
    slope = (PowerSmooth(index_peak-BWside)-PowerSmooth(index_peak+BWside))/(lambda(index_peak-BWside)-lambda(index_peak+BWside));
    Power_level_loc = PowerdB(index_peak-BWside:index_peak+BWside)-slope*(lambda(index_peak-BWside:index_peak+BWside)-lambda(index_peak+BWside))-max(PowerdB(index_peak-BWside:index_peak+BWside));
    lambda_level_loc = lambda(index_peak-BWside:index_peak+BWside);
    Power_level = PowerdB-slope*(lambda-lambda(index_peak+BWside))-max(PowerdB(index_peak-BWside:index_peak+BWside));
    else    %%% otherwise, use same fitting bandwidth from index end on.
    slope = (PowerSmooth(end-BWside*2)-PowerSmooth(end))/(lambda(end-BWside*2)-lambda(end));
    Power_level_loc = PowerdB(end-BWside*2:end)-slope*(lambda(end-BWside*2:end)-lambda(end))-max(PowerdB(end-BWside*2:end));
    lambda_level_loc = lambda(end-BWside*2:end);
    Power_level = PowerdB-slope*(lambda-lambda(end))-max(PowerdB(end-BWside*2:end));
    end
    
    %%% Resonance peak fitting
    %%% Determine "approx" resonance wavelength; as minimum
    lambda_res_guess = sum(lambda_level_loc.*((Power_level_loc-min(Power_level_loc)).^2<1e-12));
    lambda_res_index = sum((1:length(lambda)).*((lambda-lambda_res_guess).^2<1e-12));
    %%% fit around "approx" resonance with Lorentzian distribution
    if (length(lambda)-lambda_res_index)>BWside  %%% check if the peak is not "too right"
    Power_fit_dB = Power_level(lambda_res_index-BWside:lambda_res_index+BWside);
    Power_fit = 10.^(Power_fit_dB/10)*1;
    lambda_fit = lambda(lambda_res_index-BWside:lambda_res_index+BWside);
    else    %%% otherwise, use same fitting bandwidth from index end on.
    Power_fit_dB = Power_level(end-BWside*2:end);
    Power_fit = 10.^(Power_fit_dB/10)*1;
    lambda_fit = lambda(end-BWside*2:end); 
    end
    %%% lorentzian: Y(X) = P1./((X - P2).^2 + P3) + C.
    P01 = -0.001;                % init: -1
    P02 = lambda_res_guess;     
    P03 = 0.001;                % init: 0.04
    C0 = 1;                     % init: 40
    P_initial = [P01 P02 P03 C0];
    [fitting(:,peak), param(:,peak), resnorm(:,peak)] = lorentzfit(lambda_fit,Power_fit,P_initial);
    extinction_max(peak) = max(smooth(Power_fit_dB,span));
    extinction_min(peak) = min(Power_fit_dB);
    extinction_ER(peak) = max(smooth(Power_fit_dB,span))-min(Power_fit_dB);

    Power_fit_store(peak,:) = Power_fit;
    lambda_fit_store(peak,:) = lambda_fit;
end

%%% extract single-peak data; array of length #peaks
%%% extract resonance wavelength
extract.lambda_res = param(2,:);  % [nm]
%%% extract FWHM
extract.FWHM = 2*sqrt(param(3,:)); % [nm]
%%% extract extinction ratio; 
extract.extinction = extinction_ER; % [dB]
%%% extract multi-peak data; array of length #peaks-1
%%% extract FSR
extract.FSR = diff(param(2,:)); % [nm]
end