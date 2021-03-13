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
data = readtable("\Group08.xlsx", opts, "UseExcel", false); %file location on computer

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
rows = 1:88; %number of rows

dissolution_fraction = dissolution/100; %calculating dissolution fraction

% a)
figure()
plot(rows, dissolution_fraction)
title('Overview of data')
xlabel('Number of measurements')
ylabel('Dissolution_{fraction} [-]')


figure()
plot(dissolution_fraction,'o'); title('Index plot of dissolution fraction'); xlabel('Index'); ylabel('Result'); grid on; 

figure()
histogram(dissolution_fraction); title('Histogram of dissolution fraction'), xlabel('Result'), ylabel('Frequency'); grid on; 

figure()
ksdensity(dissolution_fraction); title('Density trace of dissolution fraction'), xlabel('Result'), ylabel('Density'); grid on; 

figure()
boxplot(dissolution_fraction); title('Boxplot of dissolution fraction'), xlabel('Result'); grid on; 


x = dissolution_fraction;
EDA_x = [length(x),min(x),max(x),range(x)]
EDA_y=[mean(x),median(x),length(x)] 
EDA_z=[std(x),var(x),iqr(x)] 


% b)
mean_value = mean(dissolution_fraction)
standard_diviation = std(dissolution_fraction) %standard_diviation
length_value = length(dissolution_fraction) 


x = mean(dissolution_fraction) %mean value
u = (1+standard_diviation)/(sqrt(88)) %So a 68% confidence interval leads to x +/- u, (sqrt(88) because of 88 total values and 1 because that is the z-value for a 68% confidence interval)


se = std(dissolution_fraction)/sqrt(length(dissolution_fraction)) % 68% confidence interval 



% c)
x = dissolution_fraction;
Y = log10((x/100)/(1-x/100)); %Dissolution log-ratio

% std/var(x)*(derrivative of Y with respect to x)^2
calc = std(x).*((100)/(log(10).*x.*(-x+100))).^2;

mean_value_dissolution_log = mean(calc)
se_dissolution_log = std(calc)/sqrt(length(calc))

