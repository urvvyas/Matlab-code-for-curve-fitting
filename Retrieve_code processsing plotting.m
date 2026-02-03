% Retrieves the data and plots selected systems in dB
clear all
close all

load('W15_C14_s3d1.mat','lambda','data_array', 'power_sweep')
for i=1:1:8
    Pow(:,i)=data_array{1,i}';
end
% Plot in dB - here 699.64e-6 is the input power of the laser
Plas=power_sweep.*1e-3;         % Power from laser source across the full spectra
plot(lambda, 10.*log10(Pow./Plas));
xlabel('Wavelength (nm)'); ylabel('Insertion loss (dB)')
figure;
% Plot in a normalized way for better comparison
plot(lambda, 10.*log10(Pow./max(Pow)));
xlabel('Wavelength (nm)'); ylabel('Norm. Transmission (dB)')
% 
% figure;
% for i=9:1:16
%     Pow(:,i)=data_array{1,i}';
% end
% plot(lambda, 10.*log10(Pow./max(Pow)));
% 
% figure;
% for i=17:1:24
%     Pow(:,i)=data_array{1,i}';
% end
% plot(lambda, 10.*log10(Pow./max(Pow)));
% 
% figure;
% for i=25:1:32
%     Pow(:,i)=data_array{1,i}';
% end
% plot(lambda, 10.*log10(Pow./max(Pow)));


%%%% Data in XY
% n=0;
% for i=1:1:96
% n=n+1;
% Data(:,n) = data_array{:,i};
% end