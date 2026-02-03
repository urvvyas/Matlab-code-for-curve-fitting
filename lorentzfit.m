function [fit,param,wave_fit,int_fit,P_initial] = Ring_lorentzfit(wavelength,intensity,varargin)
%%% Fit the ring spectrum to determine the full-width half-maximum
%%% bandwidth and the resonance wavelength accurately.
%%% fitting a ring spectrum using: lorentzfit(x,y)
%%% this function performs a fit according to:
%%% fit(x) = P1./((x - P2).^2 + P3) + C
%%% the function can be used as: [fit(x), param] = lorentzfit(x,y)
%%% where fit(x) are the fitted y-values
%%% param = [P1 P2 P3 C],
%%% peak amplitude = P1/P3
%%% resonance wavelength = P2
%%% FWHM ~ 2*sqrt(P3)	||||    FWHM = sqrt(8*P1*P3/(P1-2*P3*C)-4*P3)
%%% C = noise floor OR extinction coefficient; y-intercept

%%% Optional: Use a fitting bandwidth for noisy data
if nargin == 3
	BW = varargin{1};
elseif nargin > 3
    error('Too many input arguments')
else
    BW = 1;
end

[n,m] = size(wavelength);
for i=1:m
    %%% determine the (average) wavelength stepsize
    d_wave(i) = (max(wavelength(:,i))-min(wavelength(:,i)))./(n-1);
    for j=n:-1:1
        if intensity(j,i)==min(intensity(:,i))
            %%% search for the wavelength peak (minimum) at the highest 
            %%% wavelength, and determine the lower and upper wavelength 
            %%% fit bandwidth index, as well as the peak index.
            int_peak(i) = j;
            int_low_BW(i) = j-BW/2/d_wave(i);
            int_upp_BW(i) = j+BW/2/d_wave(i);
            break
        end
    end
    %%% Another option as used now: smooth outside of the peak fit BW
%     smooth_1 = smooth(intensity(1:low_BW(i),i));
%     smooth_2 = smooth(intensity(upp_BW(i):end,i));
%     intensity_smooth(:,i) = [smooth_1' intensity(low_BW(i)+1:upp_BW(i)-1,i)' smooth_2'];
%     [fit(:,i), param(:,i)] = lorentzfit(wavelength(:,i),intensity_smooth(:,i));

    %%% Determine starting values for the fit, to avoid local minima
    P_initial(:,i) = [-1 wavelength(int_peak(i)) 0.04 40];
    %%% Store the to-fit wavelength and intensity, as well as the fit
    wave_fit(:,i) = wavelength(int_low_BW(i):int_upp_BW(i),i);
    int_fit(:,i) = intensity(int_low_BW(i):int_upp_BW(i),i);
    [fit(:,i), param(:,i)] = lorentzfit(wave_fit(:,i),int_fit(:,i),P_initial(:,i));
end