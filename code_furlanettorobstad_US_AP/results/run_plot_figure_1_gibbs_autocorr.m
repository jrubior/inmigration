% this file plots the figures
clear all;
currdir = pwd;
restoredefaultpath;
% replicate Figure 1 in KM (2012)

fixed_rf=0;
label_R = 'cmy';
prior_type = 'minnesota';
ndraws  =1e6;
prior_only=0;
spec= ['rgibbs_',label_R,'prior_only_',num2str(prior_only),'prior_',prior_type,'fixed_rf_',num2str(fixed_rf),'_ndraws',num2str(ndraws)];


%L =load('matfiles/rgibbs_cmyprior_only_0prior_minnesotafixed_rf_0_ndraws100.mat')

L =load(['matfiles/',filesep,spec,'.mat']);




nhorizons = 20;%size(L.L,1);
nvar      = size(L.L,2);
nshocks   = 8;



L_q50=nan(nhorizons,nvar,nshocks); % store conjF quantile 50th
L_q16=nan(nhorizons,nvar,nshocks); % store conjF quantile 16th
L_q84=nan(nhorizons,nvar,nshocks); % store conjF quantile 84th





% figure settings
ftsizeaxis   = 11;
ftsizexlabel = 11;
ftsizetitle  = 11;
ftsizelegend = 11;
ftlinewidth  = 1.0;
medianwidth  = 1.0;


plot_uv_scale = 1;
horizon       = nhorizons-1;

if ~exist('pngfiles', 'dir')
    mkdir('pngfiles');
    disp('Created directory: pngfiles');
else
    disp('Directory pngfiles already exists');
end


close all

hFig = figure('name','autocorr');
set(hFig, 'Position', [20 20 900 900])

% Adding a main title at the top of the figure
sgtitle('Autocorrelation Functions: Variables 1 to 8', 'FontSize', 16, 'FontWeight', 'bold')

nvar1=8;
for i=1:nvar1
    for j=1:nshocks
        subplot(8,8,i+nvar1*(j-1))
        
        % Create autocorrelation plot without legend
        [acf, lags, bounds] = autocorr(squeeze(L.L(1,i,j,:)));
        
        % Manually plot the autocorrelation without using the autocorr function directly
        stem(lags, acf, 'filled', 'MarkerSize', 3)
        hold on
        plot(lags, bounds(1)*ones(size(lags)), 'r--')
        plot(lags, bounds(2)*ones(size(lags)), 'r--')
        hold off
        
        % Set axis limits similar to autocorr function
        xlim([0 lags(end)])
        ylim([-1 1])
        
        % Optional: Add grid if desired
        grid on
    end
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'L0_auto_corr_var_1to8.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'L0_auto_corr_var_1to8.png'],'-dpng');


% Adding a main title at the top of the figure
sgtitle('Autocorrelation Functions: Variables 9 to 16', 'FontSize', 16, 'FontWeight', 'bold')

nvar1=8;
for i=1:nvar1
    for j=1:nshocks
        subplot(8,8,i+nvar1*(j-1))
        
        % Create autocorrelation plot without legend
        [acf, lags, bounds] = autocorr(squeeze(L.L(1,i+8,j,:)));
        
        % Manually plot the autocorrelation without using the autocorr function directly
        stem(lags, acf, 'filled', 'MarkerSize', 3)
        hold on
        plot(lags, bounds(1)*ones(size(lags)), 'r--')
        plot(lags, bounds(2)*ones(size(lags)), 'r--')
        hold off
        
        % Set axis limits similar to autocorr function
        xlim([0 lags(end)])
        ylim([-1 1])
        
        % Optional: Add grid if desired
        grid on
    end
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'LargeSVAR_L0_auto_corr_var_9to16.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'LargeSVAR_L0_auto_corr_var_9to16.png'],'-dpng');


% Adding a main title at the top of the figure
sgtitle('Autocorrelation Functions: Variables  to 17 to 24', 'FontSize', 16, 'FontWeight', 'bold')

nvar1=8;
for i=1:nvar1
    for j=1:nshocks
        subplot(8,8,i+nvar1*(j-1))
        
        % Create autocorrelation plot without legend
        [acf, lags, bounds] = autocorr(squeeze(L.L(1,i+8*2,j,:)));
        
        % Manually plot the autocorrelation without using the autocorr function directly
        stem(lags, acf, 'filled', 'MarkerSize', 3)
        hold on
        plot(lags, bounds(1)*ones(size(lags)), 'r--')
        plot(lags, bounds(2)*ones(size(lags)), 'r--')
        hold off
        
        % Set axis limits similar to autocorr function
        xlim([0 lags(end)])
        ylim([-1 1])
        
        % Optional: Add grid if desired
        grid on
    end
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'LargeSVAR_L0_auto_corr_var_17to24.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'LargeSVAR_L0_auto_corr_var_17to24.png'],'-dpng');




% Adding a main title at the top of the figure
sgtitle('Autocorrelation Functions: Variables  to 25 to 35', 'FontSize', 16, 'FontWeight', 'bold')

nvar1=10;
for i=1:nvar1
    for j=1:nshocks
        subplot(10,8,i+nvar1*(j-1))
        
        % Create autocorrelation plot without legend
        [acf, lags, bounds] = autocorr(squeeze(L.L(1,i+8*3,j,:)));
        
        % Manually plot the autocorrelation without using the autocorr function directly
        stem(lags, acf, 'filled', 'MarkerSize', 3)
        hold on
        plot(lags, bounds(1)*ones(size(lags)), 'r--')
        plot(lags, bounds(2)*ones(size(lags)), 'r--')
        hold off
        
        % Set axis limits similar to autocorr function
        xlim([0 lags(end)])
        ylim([-1 1])
        
        % Optional: Add grid if desired
        grid on
    end
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'LargeSVAR_L0_auto_corr_var_25to35.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'LargeSVAR_L0_auto_corr_var_25to35.png'],'-dpng');




%% Responses at period 4


close all

hFig = figure('name','autocorr');
set(hFig, 'Position', [20 20 900 900])

% Adding a main title at the top of the figure
sgtitle('Autocorrelation Functions: Variables 1 to 8', 'FontSize', 16, 'FontWeight', 'bold')

nvar1=8;
for i=1:nvar1
    for j=1:nshocks
        subplot(8,8,i+nvar1*(j-1))
        
        % Create autocorrelation plot without legend
        [acf, lags, bounds] = autocorr(squeeze(L.L(4,i,j,:)));
        
        % Manually plot the autocorrelation without using the autocorr function directly
        stem(lags, acf, 'filled', 'MarkerSize', 3)
        hold on
        plot(lags, bounds(1)*ones(size(lags)), 'r--')
        plot(lags, bounds(2)*ones(size(lags)), 'r--')
        hold off
        
        % Set axis limits similar to autocorr function
        xlim([0 lags(end)])
        ylim([-1 1])
        
        % Optional: Add grid if desired
        grid on
    end
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'LargeSVAR_L4_auto_corr_var_1to8.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'LargeSVAR_L4_auto_corr_var_1to8.png'],'-dpng');


% Adding a main title at the top of the figure
sgtitle('Autocorrelation Functions: Variables 9 to 16', 'FontSize', 16, 'FontWeight', 'bold')

nvar1=8;
for i=1:nvar1
    for j=1:nshocks
        subplot(8,8,i+nvar1*(j-1))
        
        % Create autocorrelation plot without legend
        [acf, lags, bounds] = autocorr(squeeze(L.L(4,i+8,j,:)));
        
        % Manually plot the autocorrelation without using the autocorr function directly
        stem(lags, acf, 'filled', 'MarkerSize', 3)
        hold on
        plot(lags, bounds(1)*ones(size(lags)), 'r--')
        plot(lags, bounds(2)*ones(size(lags)), 'r--')
        hold off
        
        % Set axis limits similar to autocorr function
        xlim([0 lags(end)])
        ylim([-1 1])
        
        % Optional: Add grid if desired
        grid on
    end
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'LargeSVAR_L4_auto_corr_var_9to16.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'LargeSVAR_L4_auto_corr_var_9to16.png'],'-dpng');


% Adding a main title at the top of the figure
sgtitle('Autocorrelation Functions: Variables  to 17 to 24', 'FontSize', 16, 'FontWeight', 'bold')

nvar1=8;
for i=1:nvar1
    for j=1:nshocks
        subplot(8,8,i+nvar1*(j-1))
        
        % Create autocorrelation plot without legend
        [acf, lags, bounds] = autocorr(squeeze(L.L(4,i+8*2,j,:)));
        
        % Manually plot the autocorrelation without using the autocorr function directly
        stem(lags, acf, 'filled', 'MarkerSize', 3)
        hold on
        plot(lags, bounds(1)*ones(size(lags)), 'r--')
        plot(lags, bounds(2)*ones(size(lags)), 'r--')
        hold off
        
        % Set axis limits similar to autocorr function
        xlim([0 lags(end)])
        ylim([-1 1])
        
        % Optional: Add grid if desired
        grid on
    end
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'LargeSVAR_L4_auto_corr_var_17to24.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'LargeSVAR_L4_auto_corr_var_17to24.png'],'-dpng');




% Adding a main title at the top of the figure
sgtitle('Autocorrelation Functions: Variables  to 25 to 35', 'FontSize', 16, 'FontWeight', 'bold')

nvar1=10;
for i=1:nvar1
    for j=1:nshocks
        subplot(10,8,i+nvar1*(j-1))
        
        % Create autocorrelation plot without legend
        [acf, lags, bounds] = autocorr(squeeze(L.L(4,i+8*3,j,:)));
        
        % Manually plot the autocorrelation without using the autocorr function directly
        stem(lags, acf, 'filled', 'MarkerSize', 3)
        hold on
        plot(lags, bounds(1)*ones(size(lags)), 'r--')
        plot(lags, bounds(2)*ones(size(lags)), 'r--')
        hold off
        
        % Set axis limits similar to autocorr function
        xlim([0 lags(end)])
        ylim([-1 1])
        
        % Optional: Add grid if desired
        grid on
    end
end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'LargeSVAR_L4_auto_corr_var_25to35.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'LargeSVAR_L4_auto_corr_var_25to35.png'],'-dpng');





