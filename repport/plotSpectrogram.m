%%
%% Creator John Merryman - sarmap
%% Date: 02 Dec 2016
%%
%% Plot data spectrogram (time on x axis, freq. on y axis)
%%
%% usage: h = plotSpectrogram(data, plotTitle, plotXlabel, plotYlabel, ...
%%                            caxisMinMax, numSegments, freqDuration, timeDuration)
%%                          
%% data         : (input)  Data in the time domain (negative times first)
%% plotTitle, plotXlabel, plotYlabel: label strings
%% numSegments  : number of segments for short FFT calculation
%% freqDuration : Frequency duration (used for labelling)
%% timeDuration : used for labelling
%%
%% Example:
%%
%% caxisMinMax = [-80 0];
%% numSegments = 40;
%% plotTitle = 'Raw signal spectrogram';
%% plotXlabel = 'Burst time [s]';
%% plotYlabel = 'Frequency [cycles]';
%% plotSpectrogram(burstData, plotTitle, plotXlabel, plotYlabel, caxisMinMax, ...
%%    numSegments, dataNormFreqDuration, dataTimeDuration);
%%
function h = plotSpectrogram(data, plotTitle, plotXlabel, plotYlabel, caxisMinMax, numSegments, freqDuration, timeDuration)

nanMap = isnan(data);
data(nanMap) = 0;

[Sp,F,T,Ps] = spectrogram(data, numSegments);
D = flipud(fftshift(Ps,1));
%D = flipud(fftshift(Sp,1));
nli = size(D,1);
nsa = size(D,2);

figure; 
%h = imagesc(10*log10(D), caxisMinMax);
h = imagesc(10*log10(abs(D)), caxisMinMax);
colormap(winter)
c = colorbar('FontSize', 18);
c.Label.String = '[dB]';
grid on
title(plotTitle)
set(gca,'TickDir','out');

% y axis annotations                  
ylabel(plotYlabel, 'FontSize', 18)
set(gca,'YTick',[1 nli/10 2*nli/10 3*nli/10 4*nli/10 5*nli/10 6*nli/10 7*nli/10 8*nli/10 9*nli/10 nli])
%set(gca,'YTickLabel',{num2str(0.5*Ns),num2str(0.4*Ns),num2str(0.3*Ns),num2str(0.2*Ns),num2str(0.1*Ns),'0.0',num2str(-0.1*Ns),num2str(-0.2*Ns),num2str(-0.3*Ns),num2str(-0.4*Ns),num2str(-0.5*Ns)});
ytick_1  = 5*freqDuration / 10;
ytick_2  = 4*freqDuration / 10;
ytick_3  = 3*freqDuration / 10;
ytick_4  = 2*freqDuration / 10;
ytick_5  = 1*freqDuration / 10;
ytick_6  = 0;
ytick_7  = -1*freqDuration / 10;
ytick_8  = -2*freqDuration / 10;
ytick_9  = -3*freqDuration / 10;
ytick_10 = -4*freqDuration / 10;
ytick_11 = -5*freqDuration / 10;
set(gca,'YTickLabel',{num2str(ytick_1,2),num2str(ytick_2,2),num2str(ytick_3,2),num2str(ytick_4,2),...
                      num2str(ytick_5,2),num2str(ytick_6,2),num2str(ytick_7,2),num2str(ytick_8,2),...
                      num2str(ytick_9,2),num2str(ytick_10,2),num2str(ytick_11,2)});

% x axis annotations                  
xlabel(plotXlabel, 'FontSize', 18)
set(gca,'XTick',[1 nsa/10 2*nsa/10 3*nsa/10 4*nsa/10 5*nsa/10 6*nsa/10 7*nsa/10 8*nsa/10 9*nsa/10 nsa])
xtick_1  = -5*timeDuration / 10;
xtick_2  = -4*timeDuration / 10;
xtick_3  = -3*timeDuration / 10;
xtick_4  = -2*timeDuration / 10;
xtick_5  = -1*timeDuration / 10;
xtick_6  = 0;
xtick_7  = 1*timeDuration / 10;
xtick_8  = 2*timeDuration / 10;
xtick_9  = 3*timeDuration / 10;
xtick_10 = 4*timeDuration / 10;
xtick_11 = 5*timeDuration / 10;
set(gca,'XTickLabel',{num2str(xtick_1,2),num2str(xtick_2,2),num2str(xtick_3,2),num2str(xtick_4,2),...
                      num2str(xtick_5,2),num2str(xtick_6,2),num2str(xtick_7,2),num2str(xtick_8,2),...
                      num2str(xtick_9,2),num2str(xtick_10,2),num2str(xtick_11,2)}, 'FontSize', 16);