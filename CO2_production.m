clear all;close all;clc;
warning off
%% To make sure that matlab will find the functions. You must change it to your situation 
relativepath_to_generalfolder='General'; % relative reference to General folder (assumes the folder is in you working folder)
addpath(relativepath_to_generalfolder); 
%% Load Nasadatabase
TdataBase=fullfile('General','NasaThermalDatabase');
load(TdataBase);
%% Nasa polynomials are loaded and globals are set. 
%% values should not be changed. Ready for use
global Runiv Pref
Runiv=8.314472;
Pamb=1.01235e5; % Reference pressure, 1 atm!
Tamb=298.15;    % Reference Temperature
%% Some convenient units
kJ=1e3;kmol=1e3;dm=0.1;bara=1e5;kPa = 1000;kN=1000;kg=1;s=1;
%% Given conditions. For the final assignment take the ones from the specific case you are supposed to do.                           
cFuel='Gasoline';   
%% Select all species
iSp = myfind({Sp.Name},{cFuel,'O2','CO2','H2O','N2'});                      % Find indexes of these species
SpS=Sp(iSp);                                                                % Subselection of the database in the order according to {'Gasoline','O2','CO2','H2O','N2'}
%% ---------------Code------------
for percentage = [0 5 10 15]
    for load = ["no" "half" "full"]
        
        % mol_gas = (100-percentage)*volume*density/molaire massa
        mol_gas = (100-percentage)/100*4*0.75/(SpS(1).Mass*1000);                     % Density gasoline is from 0.71 to 0.77 g/ml

        % mol_ethanol = percentage*volume*density/molaire massa
        mol_ethanol = percentage/100*4*0.78945/46;                                  % Density Ethanol    0.78945 g/ml (at 20 ?C)

        x_gas = mol_gas*SpS(1).Elcomp(3); 
        y_gas = mol_gas*SpS(1).Elcomp(2); 
        z_gas = mol_gas*SpS(1).Elcomp(1); 

        x_eth = mol_ethanol*2;
        y_eth = mol_ethanol*6;
        z_eth = mol_ethanol*1;

        %% Molair masses
        C = Sp(myfind({Sp.Name},{'C'})).Mass * 1000;      % Carbon    [g/mol]
        H = Sp(myfind({Sp.Name},{'H'})).Mass * 1000;      % Hydrogen  [g/mol]
        O = Sp(myfind({Sp.Name},{'O'})).Mass * 1000;      % Oxygen    [g/mol]
        N = Sp(myfind({Sp.Name},{'N'})).Mass * 1000;      % Nitrogen  [g/mol]
        CO2 = Sp(myfind({Sp.Name},{'CO2'})).Mass * 1000;  % CO2  [g/mol]
        %% AF Gasoline blend
        M_gas = C*x_gas + H*y_gas + O*z_gas;                        % [g]
        M_eth = C*x_eth + H*y_eth + O*z_eth;                        % [g]
        M_total = M_gas + M_eth;                                    % [g] Mass fuel blend [g] for stoichiometric mix of 4 ml fuel blend
        M_CO2 = (x_gas+x_eth)*CO2;                                  % [g] Mass CO2 produced [g] for 4ml fuel blend burned
        %%
        fuel_rate_general = [14.96 11.73 8.20; 17.95 12.59 9.67; 15.71 12.95 10.82; 19.27 14.74 11.14]; % rows represent fuel rates (E0 upper row), columns represent loads (NL first column)
                                                                                             % values derived from average pressure plots experimental data                                                                                   
        if percentage == 0
            E_row = 1;
        elseif percentage == 5
            E_row = 2;
        elseif percentage == 10
            E_row = 3;
        elseif percentage == 15
            E_row = 4;
        end

        if load == "no"
            load_column = 1;
        elseif load == "half"
            load_column = 2;
        elseif load == "full"
            load_column = 3;
        end

        fuel_rate = fuel_rate_general(E_row,load_column);     % fuel rate for each blend for a specific load [4ml/s]

        %%
        %1000ml*fuel_rate(4ml/s)/4ml 
        Liter_burn_time = 1000*fuel_rate/4; %Time it takes to combust 1L fuel blend [s]
        CO2_produced = (M_CO2/4*1000)/1000; % [kg CO2/L]

        %1 liter of gasoline should produce around 2.3 kg CO2
        disp(sprintf('The motor will run for %.0f min on 1L fuel blend of E%.0f and under %s load and will produce %.3f kg CO2.', [Liter_burn_time/60, percentage, load, CO2_produced]))

        %%
        %Time per cycle
        dt_no = 0.0378;%Time per cycle, no load
        dt_half = 0.0389344545;%Time per cycle, half load
        dt_full = 0.040065;%Time per cycle, full load

        %% CO2 per cycle and CO2/hour produced
        if load == "no"
            CO2_per_cycle = M_CO2/ (fuel_rate / dt_no); % CO2 produced for 4 ml / number of cycles performed with 4 ml fuel blend  
        elseif load == "half"
            CO2_per_cycle = M_CO2/ (fuel_rate / dt_half); % CO2 produced for 4 ml / number of cycles performed with 4 ml fuel blend 
        elseif load == "full"
            CO2_per_cycle = M_CO2/ (fuel_rate / dt_full); % CO2 produced for 4 ml / number of cycles performed with 4 ml fuel blend 
        end


        Fuel_used_per_hour = (3600 / fuel_rate)*4;   % fuel used in 1 hour [ml]
        CO2_per_hour = ((Fuel_used_per_hour / 4) * M_CO2)/1000; %CO2 produced per hour [kg]
        disp(sprintf('This configuration uses %.0f ml fuel per hour and produces %.3f kg CO2 during this period.', [Fuel_used_per_hour, CO2_per_hour]))
        disp(" ")
    end
    disp(" ")
end

%% Plotting results
percentage = [0 5 10 15]; %Percentage ethanol
load = ["no load" "half load" "full load"]; % Different loads

%1 liter
Liter_burn_time = [62 49 34; 75 52 40; 65 54 45; 80 61 46] % [min]
CO2_produced = [2.407 2.407 2.407; 2.362 2.362 2.362; 2.318 2.318 2.318; 2.273 2.273 2.273] % [kg]

figure()
bar(percentage, Liter_burn_time, 'stacked')
title("Run time of 1L fuel blend")
legend(load, 'Location','northwest')
xlabel("Ethanol percentage")
ylabel('Minutes') 

figure()
bar(percentage, CO2_produced, 'stacked')
title("CO_2 produced for 1L fuel blend")
legend(load, 'Location','southeast')
xlabel("Ethanol percentage")
ylabel('Mass [kg]') 


%Per hour
Fuel_used_per_hour = [963 1228 1756; 802 1144 1489; 917 1112 1331; 747 977 1293] % [ml]
CO2_per_hour = [2.317 2.955 4.227; 1.895 2.702 3.518; 2.124 2.577 3.084; 1.698 2.220 2.938] % [kg]

figure()
bar(percentage, Fuel_used_per_hour/1000, 'stacked')
title("Fuel used per hour")
legend(load, 'Location','northeast')
xlabel("Ethanol percentage")
ylabel('Volume [L]') 

figure()
bar(percentage, CO2_per_hour, 'stacked')
title("CO_2 produced per hour")
legend(load, 'Location','northeast')
xlabel("Ethanol percentage")
ylabel('Mass [kg]') 
