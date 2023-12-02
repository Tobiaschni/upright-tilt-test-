%% Clear workspace
clear; close all; clc;


%%%%%Stand subject%%%%%%%

%%loading 
load('Sbj_5_Stand.mat','ECG','t_AUTO','RESP','Fs_Auto')

Samp_freq_base = length(ECG)

%%%preprocessing of the data 

figure
plot(t_AUTO,ECG)

%compare methods
figure
subplot(211)
plot(t_AUTO,ECG),title('ECG of the seated patient 5'),ylabel('Amplitude (mv)'),xlabel('time (s)')
subplot(212)
plot(t_AUTO,RESP),title('Respiration of the seated patient 5'),ylabel('a.u.c'),xlabel('time (s)')

% removing the mean
ECG = ECG - mean(ECG);


%removing low frequencies based drift
filter_order = 4;
fcut = 0.15;% cut-off frequency [Hz]
wc = fcut/(Fs_Auto/2); % NORMALIZED cut-off frequency
[b1,a1] = butter(filter_order,wc,'high');
figure; freqz(b1,a1,1024,Fs_Auto);

%use the filter
ECG_hp = filtfilt(b1,a1,ECG);

%%% no need to filter high frequency as the signals was already filtered tor theses frequency 

%% find R peaks

%task 7

%ECG = ECG_filtered;
[~,Rpeaks,~]=pan_tompkin(ECG,Fs_Auto,1);
RR_series = ((diff(Rpeaks))/Fs_Auto)*1000; 
t_RR = t_AUTO(Rpeaks); 
t_RR(1) = [];
RR_series = detrend(RR_series);
%removing the transitional part
n  = length(RR_series);
RR_series = RR_series(30:n);
t_RR = t_RR(30:n);

%Respiratory data
RESP = detrend(RESP);
RESPIROGRAMME = RESP(Rpeaks);
RESPIROGRAMME(1) = [];
%removing the transitional part
RESPIROGRAMME = RESPIROGRAMME(30:n);

%interpolation and resampling
n = length(t_RR);
xx = t_RR(1):0.25:t_RR(n);
RR_inter = spline(t_RR,RR_series,xx);
RESPIROGRAMME_inter = spline(t_RR,RESPIROGRAMME,xx);

%Performing the PSD with the Welch Method
F_s = 4;
N = 512;
noverlap = [];
[PSDw_ECG_stand,f3_stand] = pwelch(RR_inter,hamming(N),noverlap,[],4); 
[PSDw_RESPIRATORY_stand,f4_stand] = pwelch(RESPIROGRAMME_inter,hamming(N),noverlap,[],4); 

%%%%%%%Sit Subject%%%%%%%%%
%%loading 
load('Sbj_ 5_Sit.mat','ECG','t_AUTO','RESP','Fs_Auto')

%%%preprocessing of the data 

% removing the mean
ECG = ECG - mean(ECG);


%removing low frequencies based drift
filter_order = 4;
fcut = 0.15;       % cut-off frequency [Hz]
wc = fcut/(Fs_Auto/2); % NORMALIZED cut-off frequency
[b1,a1] = butter(filter_order,wc,'high');
figure; freqz(b1,a1,1024,Fs_Auto);



%use the filter
ECG_hp = filtfilt(b1,a1,ECG);

%%% no need to filter high frequency as the signals was already filtered tor theses frequency 

%% find R peaks

%task 7

%ECG = ECG_filtered;
[~,Rpeaks,~]=pan_tompkin(ECG,Fs_Auto,1);
RR_series = ((diff(Rpeaks))/Fs_Auto)*1000; 
t_RR = t_AUTO(Rpeaks); 
t_RR(1) = [];
RR_series = detrend(RR_series);
%removing the transitional part
n  = length(RR_series);
RR_series = RR_series(30:n);
t_RR = t_RR(30:n);

%Respiratory data
RESP = detrend(RESP);
RESPIROGRAMME = RESP(Rpeaks);
RESPIROGRAMME(1) = [];
%removing the transitional part
RESPIROGRAMME = RESPIROGRAMME(30:n);

% interpolation and resampling
n = length(t_RR);
xx = t_RR(1):0.25:t_RR(n);
RR_inter = spline(t_RR,RR_series,xx);
RESPIROGRAMME_inter = spline(t_RR,RESPIROGRAMME,xx);

%Performing the PSD with the Welch Method
F_s = 4;
N = 512;
noverlap = [];
[PSDw_ECG_sit,f3_sit] = pwelch(RR_inter,hamming(N),noverlap,[],4); %10 windows, 50% overlap - Welch
[PSDw_RESPIRATORY_sit,f4_sit] = pwelch(RESPIROGRAMME_inter,hamming(N),noverlap,[],4); %10 windows, 50% overlap - Welch



figure
plot(f3_stand,PSDw_ECG_stand)
xlim([0,0.5]),title('PSD of HRV'),xlabel('Hz'),ylabel('msec^2/Hz');


figure
plot(f4_stand,PSDw_RESPIRATORY_stand)
xlim([0,0.5]),title('PSD of RESPIRATION'),xlabel('Hz'),ylabel('a.u.c^2/Hz');

%compare methods
figure
subplot(211)
plot(f3_sit,PSDw_ECG_sit,f3_stand,PSDw_ECG_stand),title('Comparison HRV PSD Sit/Stand'),ylabel('Power/f'),xlabel('Hz'),xlim([0,0.5])
legend('Sit','Stand')
subplot(212)
plot(f4_sit,PSDw_RESPIRATORY_sit,f4_stand,PSDw_RESPIRATORY_stand),title('Comparison RESPIROGRAM PSD sit/stand'),ylabel('Power/f'),xlabel('Hz')
legend('Sit','Stand')

figure
plot(f3_sit,PSDw_ECG_sit,f3_stand,PSDw_ECG_stand),title('Comparison HRV PSD Sit/Stand'),ylabel('Power/f'),xlabel('Hz'),xlim([0,0.5])
legend('Sit','Stand')

% identification LF / HF

i_min = 1;
while f3_sit(i_min) < 0.05
    i_min = i_min+1;
end
 

i_max = 1;
while f3_sit(i_max) < 0.14
    i_max = i_max+1;
end

i_mmax = 1;
while f3_sit(i_mmax) < 0.4
    i_mmax = i_mmax+1;
end

puissance_tot_LF_sit = trapz(f3_sit(i_min:i_max),PSDw_ECG_sit(i_min:i_max));
puissance_tot_HF_sit = trapz(f3_sit(i_max:i_mmax),PSDw_ECG_sit(i_max:i_mmax));

sprintf( "LF %s",puissance_tot_LF_sit)
sprintf( "HF %s",puissance_tot_HF_sit)
sprintf( "le ratio sit est %s",puissance_tot_LF_sit/puissance_tot_HF_sit)



%cross spectrale analysis
Fs= 4;
nfft = length(RR_inter); 
x = RR_inter;
y = RESPIROGRAMME_inter;


%Cross-spectrum amplitude
[Pxy,f] = cpsd(x,y,hamming(N),256,nfft,Fs); 
figure
plot(f,abs(Pxy)),xlim([0,0.5]),title('Cross Spectral Amplitude of HRV and Respirogram'),xlabel('Hz'),ylabel('Pxy')


%Cross-spectrum phase
phase = angle(Pxy);
figure
plot(f,phase),xlim([0,0.5]),title('Cross Spectral Phase of HRV and Respirogram'),xlabel('Hz'),ylabel('phi (rad)');


% %Quadratic coherence
[Cxy,f] = mscohere(x,y,hamming(256),128,nfft,Fs);
figure
plot(f,Cxy),xlim([0,0.5]),xlabel('Hz'),ylabel('k^2');
hold on
plot([0,0.5],[0.5,0.5])
hold on; 
yyaxis right;
plot(f,phase,'--'),xlim([0,0.5]),ylim([-3.14,3.14]),title('Coherence function and phase of the cross spectral analysis')
legend('amplitude','','phase')



  