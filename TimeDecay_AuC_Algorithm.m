clear
clc
warning('off','all')
% Input directory for the folder that data is saved to
path = 'C:\YOUR DIRECTORY';    

% Input file name beginning with forward slash
file = 'CNMF OUTPUT FILE NAME';   
filenamestr = 'FIGURE TITLE STRING'; %cannot use udnerscores due to MATLAB syntax  

%


max_tolerance = 0.9;        % Logic gating for findpeaks function, can modify if incorrect peaks are being identified
fc = 0.05;                  % Cut off frequency; fc*sampling frequency = true cut off freq
fps = 1.41;                 % Parameter of recording device
low_thresh = 0.15;          % Threshold for end of time decay calculation (85% of peak)
nic = 5;                    % Enter time when you injected drug #1 in minutes (Originally Nicotine)
kcl = 20;                   % Enter time when you injected drug #2 in minutes (Originally KCl)
insig_thresh = 0.75;        % If Peak signal does exceed this threshold, not included in plotting
                            

%-----------------Changes not needed below----------------                           
filepath = [path '\' file '.mat'];
sg_combined = load(filepath);
figurestrname = 'Output_Figures';
figurestr = [figurestrname '_' file];
mkdir(path,figurestr);
folderpath = [path '\' figurestr];

% Function outputting Excel matrix
[trace_table] = Decay_AuC(sg_combined,fps, fc, max_tolerance,low_thresh,filenamestr,nic,kcl,folderpath, insig_thresh);

% Code to export trace_table as Excel file within defined directory
trace_tablestr = 'Algorithm_Outputs';
write_trace_table = [path '\' file '_' trace_tablestr '.xlsx'];
writetable(trace_table,write_trace_table)

%email balthazor@ohsu.edu or ntbalthazor@gmail.com for further questions with data processing or
%potential edits



