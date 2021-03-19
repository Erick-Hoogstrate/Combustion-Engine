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
data = readtable("C:\Users\20192651\Documents\Year 2\Q3\Combustion Engine 4GB10\AnaMeDa\Group08.xlsx", opts, "UseExcel", false); %file location on computer

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

%% Task 2.d

% scatter(blend_speed,dissolution)                    %make scatter plots to see correlation
% scatter(force,dissolution)
%scatter(spray_rate,dissolution);
% scatter(atomizer_pressure,dissolution)


a=corrcoef(dissolution,blend_speed); %0.2116          %calculating correlation coefficient
b=corrcoef(dissolution,force); %0.1597
c=corrcoef(dissolution,spray_rate); %-0.3009
d=corrcoef(dissolution,atomizer_pressure); %0.1009

Simpel_reg = fitlm(spray_rate,dissolution)          %fitting a simple regression model
figure()
plot(Simpel_reg)

cof=coefCI(Simpel_reg)                              %confirmation that B0 and B1 are significant in 95% confidence interval

%% Task 2.e residual plot
%plotResiduals(Simpel_reg, 'fitted');
%plotResiduals(Simpel_reg, 'caseorder');

%plotResiduals(Simpel_reg)
%plotResiduals(Simpel_reg,'probability')

%% Task 2.f
disspray=table(spray_rate, dissolution);    
quad_reg=fitlm(disspray,'dissolution ~ spray_rate + spray_rate^2')  %make quadratic regression model
figure()
plot(quad_reg)

%Rsquared higher so better

%% Task 3.g

 CIFcn = @(dissolution,p)prctile(x,abs([0,100]-(100-p)/2));
 % 68% interval for Dissolution Fraction (Direct approach):
 CI68 = CIFcn(dissolution,68); % 68% confidence interval

%Creating scatter plots of different variables

figure()
subplot(3, 2, 1);
plot(blend_speed,force, 'o')
xlabel('Blend speed')
ylabel('force')

subplot(3, 2, 2);
plot(blend_speed,spray_rate, 'o')
xlabel('Blend speed')
ylabel('Spray rate')

subplot(3, 2, 3);
plot(blend_speed,atomizer_pressure, 'o')
xlabel('Blend speed')
ylabel('Atomizer pressure')

subplot(3, 2, 4);
plot(force,spray_rate, 'o')
xlabel('Force')
ylabel('Spray rate')

subplot(3, 2, 5);
plot(force,atomizer_pressure, 'o')
xlabel('Force')
ylabel('Atomizer pressure')

subplot(3, 2, 6);   
plot(spray_rate,atomizer_pressure, 'o')
xlabel('Spray rate')
ylabel('Atomizer pressure')

e=corrcoef(blend_speed, force);
f=corrcoef(blend_speed, spray_rate);
g=corrcoef(blend_speed, atomizer_pressure);
h=corrcoef(force, spray_rate);
i=corrcoef(force, atomizer_pressure);
j=corrcoef(spray_rate, atomizer_pressure);

% Checking for normality with histograms and qq-plots.  
figure();
grid on;
subplot(2, 2, 1);
CI68v = CIFcn(blend_speed,68);
histfit(blend_speed);
title('Histogram of blend speed');
xline(CI68v(1));
xline(CI68v(2));

subplot(2, 2, 2);
CI68v = CIFcn(force,68);
histfit(force);
title('Histogram of force');
xline(CI68v(1));
xline(CI68v(2));

subplot(2, 2, 3);
CI68v = CIFcn(spray_rate,68);
histfit(spray_rate);
title('Histogram of spray rate');
xline(CI68v(1));
xline(CI68v(2));

subplot(2, 2, 4);
CI68v = CIFcn(atomizer_pressure,68);
histfit(atomizer_pressure);
title('Histogram of atomizer pressure');
xline(CI68v(1));
xline(CI68v(2));

figure();
grid on;
subplot(2,2,1);
qqplot(blend_speed)
title('QQ-plot of blend speed');

subplot(2,2,2);
qqplot(force)
title('QQ-plot of force');

subplot(2,2,3);
qqplot(spray_rate)
title('QQ-plot of spray rate');

subplot(2,2,4);
qqplot(atomizer_pressure)
title('QQ-plot of atomizer pressure');