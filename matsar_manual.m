%% Stochastic gradient descent via Constant Update, AdaGrad or AdaDelta.
% For debugging if layout-introduced parastics matters to performance,
% use either
% 1. large mu and smaller tend (<1e-4). The SNR will not be > 62DB (or
% >10ENOB) but it shows quickly if calibration is possible.
% 2. Constant mu (1e-6) and tend (>1e-3). The SNR may be > 62DB depending
% on the input.
%
% To obtain final weights, use either
% 1. AdaDelta and tend > 1e-2
% 2. AdaGrad and mu < 1e-3 and tend > 1e-2
% 3. Constant mu < 1e-7 and tend > 1e-3
%
% In all cases, AdaDelta guarantees best result if time permits.
%
% Capacitors are changed via the following modifier
% Cactual = C.*Cpar_mult+Cpar_add;
% where each element  Cpar_mult and Cpar_add applies each capacitor (MSB
% first).
%
% Use C_mask to disable or enable elements.
clear;

mode = 'AdaDelta';  % Enables either, Constant, AdaGrad or AdaDelta
mu = 1e-4; % Initial learning rate.

tend = 1e-3;
C_mask = logical([1,1,1,1,1,1,1,1,1,1,1,1,1]);
Cpar_mult = [1,1,1,1,1,1,1,1,1,1,1,1,1];
Cpar_add = [-2.7900e-15, -2.1190e-15, -2.7720e-15,...
    -1.3180e-15, -5.3390e-16, -2.7950e-16, ...
    -1.0410e-16, 1.6160e-17, 6.1750e-17, ...
    9.017e-17, 5.608e-17, 7.629e-17, 1.0520e-16];
%Cpar_add = [0,0,0,0,0,0,0,0,0,0,0,0,0];

Cmin = 1e-15;
FullScale = 1.2;
Carray = [4.72e-13, 2.36e-13, 1.36e-13,...
    7.6e-14, 4.4e-14, 2.4e-14,...
    1.6e-14, 8e-15, 4e-15,...
    4e-15, 2e-15, 1e-15, 1e-15];
M = sum(C_mask)-1;

C = Carray(C_mask);
CT = sum(C);
C_mis=Carray.*Cpar_mult+Cpar_add;
C_mis=C_mis(C_mask);
CT_mis = sum(C_mis);

Rdn = zeros(1,M);
for n = 1:M
    Rdn(n) = (sum(C(n+1:end))-C(n));
end
Rdn_Percent = Rdn./Carray(n);
Rdn = Rdn./Cmin;

kT = 1.38064852e-23*350; %T=350K

fs = 204.8e6; ts = 1/fs;
numsamp = floor(tend/ts);
t = linspace(0,tend-ts,numsamp);
rbw = fs/numsamp;

f = 1e6; w = 2*pi*f; Vcm = FullScale/2;
Input = (FullScale*1/2)*(sin(w*t))+Vcm;

W = ones(numsamp,M);
Regs = zeros(numsamp,M);
Regs_mis = zeros(numsamp,M);
X = zeros(numsamp,M);
err = zeros(1,numsamp);
Gamma = 0.7;
Eg = 0; EdW  = 0;
for n = 1:numsamp
    for ii = 1:M
        Regs(n,ii) = 1; Regs_mis(n,ii) = 1;
        Regs(n,ii) = int8(Input(n) > (FullScale * (C(1:M) * Regs(n,:)')/CT));
        Regs_mis(n,ii) = int8(Input(n) > (FullScale * (C_mis(1:M) * Regs_mis(n,:)')/CT_mis) );
    end
    % Run LMS on current sample
    % Update rule W[n+1] = W[n]+(Y-WX)X
    X(n,:) = Regs_mis(n,:) .* C(1:M) /Cmin;
    Y = Regs(n,:) * C(1:M)'/Cmin; %0-1023
    err(n) = Y - (W(n,:) * X(n,:)'); %scalar
    grad = err(n)*X(n,:); %vector
    if n<numsamp
        switch (mode)
            case 'AdaDelta'
                Eg = Gamma*Eg+(1-Gamma)*grad.^2;
                if n>1
                    dWs = (W(n,:)-W(n-1,:)).^2;
                    EdW = Gamma*EdW + (1-Gamma)*dWs;
                    mup = sqrt(1e-12+EdW)./sqrt(1e-12+Eg);
                else
                    mup = mu./sqrt(1e-12+Eg);
                end
            case 'AdaGrad'
                Eg = Eg + grad.^2;
                mup = mu./sqrt(1e-12+Eg);
            case 'Constant'
                mup = mu;
        end
        W(n+1,:) = W(n,:) + mup.*grad;
    end
end
% Apply final weights to calculate actual decimal readout
codes_ideal = zeros(1,numsamp);
codes_precal = zeros(1,numsamp);
codes_cal = zeros(1,numsamp);
for n = 1:numsamp
    codes_ideal(n)  = Regs(n,:) * C(1:M)'/Cmin;
    codes_cal(n) = W(numsamp,:)*(Regs_mis(n,:) .* C(1:M) /Cmin)';
    codes_precal(n) = Regs_mis(n,:) * C(1:M)'/Cmin;
end


%% Plot section
Windowing = 1;
ShowWeights = 1; % Shows weight history
ShowError = 1;  % Shows error history
ShowSE = 1;
ShowWave = 1; % Shows time domain waveforms.
ShowSpectrum = 1; % Shows spectra.
ShowSNDR = 1; % Shows SNR
ShowSFDR = 1; % Shows SFDR
ShowENOB = 1; % Shows ENOB
ShowDNL = 0;
ShowHistogram = 0;
ShowINL = 0;

dbv=@(v) 20.*log10(abs(v));
noise = 4*kT/CT/numsamp;
noiseDB = 10*log10(noise);
clc;
close all;

if ShowWave
    figure
    subplot(2,1,1);
    plot(t,codes_precal,'-'); hold on;
    plot(t,codes_ideal,':');
    xlim([1/f,3/f]); title('Precal');
    hold  off; legend('PreCal','Ideal');
    subplot(2,1,2);
    plot(t,codes_cal,'-'); hold  on;
    plot(t,codes_ideal,':');
    xlim([1/f,3/f]); title('Cal');
    hold  off; legend('Cal','Ideal');
end

if ShowError
    figure
    plot(err); title('Error');
    xlabel('Number of Samples');
    ax1 = gca;
    ax1.XGrid = 'on';
    ax1.YGrid = 'on';
    ax1.Position = ax1.Position + [0 0.15 0 -0.15];
    ax1.YLabel.String = 'LSB';
    ax1.YLim = [-10,10];
    ax1.XLim = [1,numsamp];
    ax2 = axes('position', (ax1.Position .* [1 1 1 1]) + [0 -0.15 0 0], 'color', 'none', 'YColor','none','linewidth', 1);
    ax2.XLim = [min(t),max(t)];
    ax2.XLabel.String = 'Time (s)';
end

if ShowSE
    figure
    plot(err.^2); title('Squared Error');
    xlabel('Number of Samples');
    ax1 = gca;
    ax1.XGrid = 'on';
    ax1.YGrid = 'on';
    ax1.Position = ax1.Position + [0 0.15 0 -0.15];
    ax1.YLabel.String = 'LSB';
    ax1.YLim = [0,10];
    ax1.XLim = [1,numsamp];
    ax2 = axes('position', (ax1.Position .* [1 1 1 1]) + [0 -0.15 0 0], 'color', 'none', 'YColor','none','linewidth', 1);
    ax2.XLim = [min(t),max(t)];
    ax2.XLabel.String = 'Time (s)';
end

% calculates ideal voltage level
VbinsIdeal = zeros(1,1024);
VbinsPrecal = zeros(1,1024);
VbinsCal = zeros(1,1024);
HisIdeal = zeros(1,1024);
HisPrecal = zeros(1,1024);
HisCal = zeros(1,1024);
P = zeros(1,1024);

V = @(c) (c-512)*FullScale./1024;
for ii = 1:1024
    % Doernberg
    HisIdeal(ii) = nnz(codes_ideal<ii & codes_ideal>=ii-1);
    HisPrecal(ii) = nnz(codes_precal<ii & codes_precal>=ii-1);
    HisCal(ii) = nnz(codes_cal<ii & codes_cal>=ii-1);
    %P(ii) = (asin(V(ii)*2/FullScale)-asin(V(ii-1)*2/FullScale))/pi;
end

for ii = 1:1024
    VbinsIdeal(ii) = -cos(pi*sum(HisIdeal(1:ii))/numsamp);
    VbinsPrecal(ii) = -cos(pi*sum(HisPrecal(1:ii))/numsamp);
    VbinsCal(ii) = -cos(pi*sum(HisCal(1:ii))/numsamp);
end

if ShowHistogram
    figure
    subplot(3,1,1); 
    bar(HisIdeal./numsamp,'LineWidth',2); title('Ideal');
    ylabel('Normalized'); xlabel('Code'); ylim([0 0.025]);
    ax=gca; ax.TickLength = [0.005,0.025];
    xlim([-10,1033]); grid on;
    subplot(3,1,2); 
    bar(HisPrecal./numsamp,'LineWidth',2); title('Precal');
    ylabel('Normalized'); xlabel('Code');
    ax=gca; ax.TickLength = [0.005,0.025]; ylim([0 0.025]);
    xlim([-10,1033]); grid on;
    subplot(3,1,3);
    bar(HisCal./numsamp,'LineWidth',2); title('Cal');
    ylabel('Normalized'); xlabel('Code');
    ax=gca; ax.TickLength = [0.005,0.025]; ylim([0 0.025]);
    xlim([-10,1033]); grid on;
end

if ShowDNL
    figure
    subplot(3,1,1);
    DNLideal = (diff(VbinsIdeal)/(VbinsIdeal(1023)-VbinsIdeal(1))-Cmin/CT)*512;
    plot(DNLideal);
    xlim([0,1024]); title('Ideal'); ylim([-2,2]);
    grid on;
    ylabel('DNL'); xlabel('Code');
    subplot(3,1,2);
    DNLprecal = (diff(VbinsPrecal)/(VbinsPrecal(1023)-VbinsPrecal(1))-Cmin/CT)*512;
    plot(DNLprecal);
    grid on;
    xlim([0,1024]); title('Precal'); ylim([-2,2]);
    ylabel('DNL'); xlabel('Code');
    subplot(3,1,3);
    DNLcal = (diff(VbinsCal)/(VbinsCal(1023)-VbinsCal(1))-Cmin/CT)*512;
    plot(DNLcal);
    grid on;
    xlim([0,1024]); title('Cal'); ylim([-2,2]);
    ylabel('DNL'); xlabel('Code');
end

if ShowINL
    figure
    
    subplot(3,1,1);
    INLideal = zeros(1,1023);
    for ii = 2:1024
        INLideal(ii) = sum(DNLideal(1:ii-1));
    end
    plot(INLideal);
    grid on;
    xlim([0,1024]); title('Ideal');
    ylabel('INL'); xlabel('Code');
    
    subplot(3,1,2);
    INLprecal = zeros(1,1023);
    for ii = 2:1024
        INLprecal(ii) = sum(DNLprecal(1:ii-1));
    end
    plot(INLprecal);
    grid on;
    xlim([0,1024]); title('Precal');
    ylabel('INL'); xlabel('Code');
    
    subplot(3,1,3);
    INLcal = zeros(1,1023);
    for ii = 2:1024
        INLcal(ii) = sum(DNLcal(1:ii-1));
    end
    plot(INLcal);
    grid on;
    xlim([0,1024]); title('Cal');
    ylabel('INL'); xlabel('Code');
end

Win = ones(1,numsamp);
if Windowing
    Win = blackman(numsamp)';
end
fftr = abs(fft(codes_ideal.*1.2/1024.*Win))/numsamp;
fftr = max(noise,fftr(1:floor(numsamp/2)));

if ShowSNDR
    sndrIdeal = dbv(SNDR(fftr));
    fprintf('SNDR(Ideal) = %f\n',sndrIdeal);
end
if ShowSFDR
    sfdrIdeal = dbv(SFDR(fftr));
    fprintf('SFDR(Ideal) = %f\n',sfdrIdeal);
end
if ShowENOB
    enobIdeal = (dbv(SNDR(fftr))-1.76)./6.02;
    fprintf('ENOB(Ideal) = %f\n',enobIdeal);
end
if ShowSpectrum
    figure
    f=(0:numsamp/2-1)*fs/numsamp;
    subplot(3,1,1);
    plot(f,dbv(fftr));  title('Ideal');
    ylim([-120,0]);
end

fftr = abs(fft(codes_precal.*1.2/1024.*Win,numsamp)/numsamp);
fftr = max(noise,fftr(1:floor(numsamp/2)));
if ShowSNDR
    sndrPrecal = dbv(SNDR(fftr));
    fprintf('SNDR(Precal) = %f\n',sndrPrecal);
end
if ShowSFDR
    sfdrPrecal = dbv(SFDR(fftr));
    fprintf('SFDR(Precal) = %f\n',sfdrPrecal);
end
if ShowENOB
    enobPrecal = (dbv(SNDR(fftr))-1.76)./6.02;
    fprintf('ENOB(Precal) = %f\n',enobPrecal);
end
if ShowSpectrum
    subplot(3,1,2);
    plot(f,dbv(fftr));  title('Precal');
    ylim([-120,0]);
end


fftr = abs(fft(codes_cal.*1.2/1024.*Win,numsamp)/numsamp);
fftr = max(noise,fftr(1:floor(numsamp/2)));
if ShowSNDR
    sndrCal = dbv(SNDR(fftr));
    fprintf('SNDR(Cal) = %f\n',sndrCal);
end
if ShowSFDR
    sfdrCal = dbv(SFDR(fftr));
    fprintf('SFDR(Cal) = %f\n',sfdrCal);
end
if ShowENOB
    enobCal = (dbv(SNDR(fftr))-1.76)./6.02;
    fprintf('ENOB(Cal) = %f\n',enobCal);
end
if ShowSpectrum
    subplot(3,1,3);
    plot(f,dbv(fftr));  title('Cal');
    ylim([-120,0]);
end

if ShowWeights
    for k = 1:M
        l=mod(k-1,6);
        if l == 0
            figure;
        end
        subplot(3,2,l+1); plot(t,W(:,k));
        title(['W',num2str(M-k)]);
    end
end
%lows =  Regs_mis(Regs_mis(:,1)==0,:);