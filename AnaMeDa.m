clc; clear all; close all

%% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 41);

% Specify sheet and range
opts.Sheet = "Sheet2";
opts.DataRange = "A2:AO89";

% Specify column names and types
opts.VariableNames = ["Lot", "BlendSpeed", "Compressor", "Force", "SprayRate", "AtomizerPressure", "Dissolution", "obs1", "obs2", "obs3", "obs4", "obs5", "obs6", "obs7", "obs9", "obs10", "obs11", "obs12", "obs13", "obs14", "obs15", "obs16", "obs17", "obs18", "obs19", "obs20", "obs21", "obs22", "obs23", "obs24", "obs25", "obs26", "obs27", "obs28", "obs29", "obs30", "obs31", "obs32", "obs33", "obs34", "obs35"];
opts.SelectedVariableNames = ["Lot", "BlendSpeed", "Compressor", "Force", "SprayRate", "AtomizerPressure", "Dissolution", "obs1", "obs2", "obs3", "obs4", "obs5", "obs6", "obs7", "obs9", "obs10", "obs11", "obs12", "obs13", "obs14", "obs15", "obs16", "obs17", "obs18", "obs19", "obs20", "obs21", "obs22", "obs23", "obs24", "obs25", "obs26", "obs27", "obs28", "obs29", "obs30", "obs31", "obs32", "obs33", "obs34", "obs35"];
opts.VariableTypes = ["string", "double", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, [1, 3], "EmptyFieldRule", "auto");

% Import the data
data = readtable("C:\Users\20192420\Documents\Q3 Y2\AnaMeDa\Group08.xlsx", opts, "UseExcel", false); %file location on computer

lot               = data.Lot;
blend_speed       = data.BlendSpeed;
compressor        = data.Compressor;
force             = data.Force;
spray_rate        = data.SprayRate;
atomizer_pressure = data.AtomizerPressure;
dissolution       = data.Dissolution;

% Clear temporary variables
clear opts

%% 
rows = 1:88 %number of rows

dissolution_fraction = dissolution/100 %calculating dissolution fraction

figure()
plot(rows, dissolution_fraction)
title('Overview of data')
xlabel('all different measurements')
ylabel('dissolution_fraction')

standard_diviation = std(dissolution_fraction) %standard_diviation
x = mean(dissolution_fraction) %mean value

u = (1+standard_diviation)/(sqrt(88)) %So a 68% confidence interval leads to x +/- u, (sqrt(88) because of 88 total values and 1 because that is the z-value for a 68% confidence interval)

