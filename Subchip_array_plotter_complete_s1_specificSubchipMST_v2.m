%%% Plot measurement results per row in groups
clear;
close all;

%%% Adding data path, if not placed in same folder:
% addpath('H:\MST_Austausch_double\Analysis Tools\Subchiparray Plotter');
load('W15_C14_s3d1.mat','lambda','data_array')
%%% subchip 1 analysis for specific devices of a specific chip
chip = 14;                          %%% select chip, this may be extended to also use the subchip and row (d) parameters 
FSR = 4.5;                          %%% approx. FSR of used devices, this may be optimzed later to be adaptive


%%% different version to call data
% devices = [1 2 4 5 6 7 8]; % Create the array of devices manually
% devices = [92 93 94 95 96]; 
devices = linspace(1,16,16); % Create the array of devices with linespace
d_index = 1;

%%% load data for chip
data=load(strcat('W15_C',num2str(chip),'_s3d1'));
close all

%%% Adjust data: 'normalize'
for d=1:length(data.data_array)
    data.data_array{d} = data.data_array{d}/max(data.data_array{d})*1e-3;
end
%%% plot device data
fig1=figure('Units','normalized','Position',[0 0 1 1]);
plot(data.lambda,10.*log10(data.data_array{devices(1)}./1e-3),'k','linewidth',2)
hold on
for jj=devices(2:end)
plot(data.lambda,10.*log10(data.data_array{jj}./1e-3),'k','linewidth',2)
end
hold off
set(gca,'FontSize',45)
set(gca,'linewidth',4)
set(gca, 'FontName', 'Times')
xlabel('Wavelength (nm)')
ylabel('Norm. Transmission (dB)')

%%% fit resonance peaks
for ii=1:length(data)
    for jj=1:length(devices)
        %%% fit peaks right to left
[Power_fit_store{jj,ii},lambda_fit_store{jj,ii},fitting{jj,ii},...
    extract{jj,ii},param{jj,ii}] = ResonancePeakAnalysis_lin_back(...
    data.data_array{devices(jj)},data.lambda,FSR);
        
        if max(extract{jj,ii}.extinction(end-1))<3
        disp(strcat('Bad ER or poor fit for (jj,ii)=(',num2str(jj),',',num2str(ii),'). Refit'))
        %%% fit peaks left to right
[Power_fit_store{jj,ii},lambda_fit_store{jj,ii},fitting{jj,ii},...
    extract{jj,ii},param{jj,ii}] = ResonancePeakAnalysis_lin(...
    data.data_array{devices(jj)},data.lambda,FSR);
        end
    end
end

jj = d_index;
fig2=figure('Units','normalized','Position',[0 0 1 1]);
plot(data.lambda,double(10.*log10(data.data_array{devices(jj)}*1e3)),'k','LineWidth',2)
hold on
for pp=1:length(lambda_fit_store{jj}(:,1))
plot(lambda_fit_store{jj}(pp,:),10*log10(Power_fit_store{jj}(pp,:)),'b','LineWidth',2,'LineStyle','none','Marker','o')
plot(lambda_fit_store{jj}(pp,:),10*log10(fitting{jj}(:,pp)),'r','LineWidth',2)
end
hold off
set(gca,'FontSize',45)
set(gca,'linewidth',4)
set(gca, 'FontName', 'Times')
legend('data','leveled data','fit','location','SouthEast')
xlabel('Wavelength (nm)')
ylabel('Norm. Transmission (dB)')
title('Lorentz fit to leveled MRR data')

fig3=figure('Units','normalized','Position',[0 0 1 1]);
plot(extract{1}.lambda_res,extract{1}.extinction,'k','LineWidth',2,'LineStyle','none','Marker','o')
hold on
for jj=2:length(devices)
plot(extract{jj}.lambda_res,extract{jj}.extinction,'k','LineWidth',2,'LineStyle','none','Marker','o')
end
hold off
% axis([1480 1560 -30 -10])
set(gca,'FontSize',45)
set(gca,'linewidth',4)
set(gca, 'FontName', 'Times')
xlabel('Resonance wavelength (nm)')
ylabel('Extinction ratio (dB)')