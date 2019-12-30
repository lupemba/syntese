close all; clear all; clc;
data1 = csvread("real_original.csv");
data2 = csvread("imag_original.csv");
data_original = data1 + 1i*data2;

data1 = csvread("real_deramped.csv");
data2 = csvread("imag_deramped.csv");
data_deramped = data1 + 1i*data2;

data1 = csvread("real_resample.csv");
data2 = csvread("imag_resample.csv");
data_resample = data1 + 1i*data2;

data1 = csvread("real_reramped.csv");
data2 = csvread("imag_reramped.csv");
data_reramped = data1 + 1i*data2;

caxisMinMax = [-80 0];
numSegments = 40;
plotTitle = 'Raw signal spectrogram';
plotXlabel = 'Burst time [s]';
plotYlabel = 'Frequency [cycles]';
freqDuration = 1;
timeDuration = 3;
%%
plotSpectrogram(data_original, plotTitle, plotXlabel, plotYlabel, caxisMinMax, numSegments, freqDuration, timeDuration)
%%
plotSpectrogram(data_deramped, plotTitle, plotXlabel, plotYlabel, caxisMinMax, numSegments, freqDuration, timeDuration)
%%
plotSpectrogram(data_resample, plotTitle, plotXlabel, plotYlabel, caxisMinMax, numSegments, freqDuration, timeDuration)
%%
plotSpectrogram(data_reramped, plotTitle, plotXlabel, plotYlabel, caxisMinMax, numSegments, freqDuration, timeDuration)