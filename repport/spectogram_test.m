close all; clear all; clc;
data1 = csvread("original_real.csv");
data2 = csvread("original_imag.csv");
data_original = data1 + 1i*data2;

data1 = csvread("deramped_real.csv");
data2 = csvread("deramped_imag.csv");
data_deramped = data1 + 1i*data2;

data1 = csvread("chirp_real.csv");
data2 = csvread("chirp_imag.csv");
data_resample = (data1 + 1i*data2);

data1 = csvread("reramped_real.csv");
data2 = csvread("reramped_imag.csv");
data_reramped = data1 + 1i*data2;

caxisMinMax = [-80 0];
numSegments = 40;
plotTitle = 'Raw signal spectrogram';
plotXlabel = 'Azimuth time, wrt. mid-burst time [s]';
plotYlabel = 'Normalised Frequency [.]';
freqDuration = 1;
timeDuration = 3;
%%
plotSpectrogram(data_original, '', plotXlabel, plotYlabel, caxisMinMax, numSegments, freqDuration, timeDuration)
%%
plotSpectrogram(data_deramped, '', plotXlabel, plotYlabel, caxisMinMax, numSegments, freqDuration, timeDuration)
%%
plotSpectrogram(data_resample, '', plotXlabel, plotYlabel, caxisMinMax, numSegments, freqDuration, timeDuration)
%%
plotSpectrogram(data_reramped, '', plotXlabel, plotYlabel, caxisMinMax, numSegments, freqDuration, timeDuration)