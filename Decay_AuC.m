function [trace_table] = Decay_AuC(sg_combined,fps,fc,max_tolerance,low_scale,filenamestr,nic,kcl, folderpath, insig_thresh)
%DECAY_AOC Summary of this function goes here
%   Detailed explanation goes here
num_neurons = length(sg_combined.F_raw(:,1));
removed = zeros(num_neurons,1);
tracelength = length(sg_combined.F_raw(1,:));


%defining nicotine and kcl injection boundary that divides peaks 
nic_bound = int32((nic-0.5)*60*fps);
kcl_bound = int32((kcl-0.5)*60*fps);

%preallocating output arrays
vertgap = cell(num_neurons,1);
nicnamecol = cell(num_neurons,1);
nicnamecol(1) = cellstr('Nicotine');
kclnamecol = cell(num_neurons,1);
kclnamecol(1) = cellstr('KCl');
neuroncol = cell(num_neurons,1);
maxes = zeros(num_neurons,1);
filtset = zeros(num_neurons,tracelength);
for d=1:num_neurons
    num = sprintf('%i',d);
    neuroncol(d) = cellstr(num);
    axis_set_raw = sg_combined.F_raw(d,:)';
    axis_set = axis_set_raw - min(axis_set_raw);
    [b,a] = butter(4,fc,'low');
    filtaxis = filtfilt(b,a,axis_set);
    filtaxis_1 = filtaxis - min(filtaxis);
    filtset(d,:) = filtaxis_1;
    maxes(d) = max(filtaxis_1);
end
nic_max_dffs = zeros(num_neurons,1);
nic_decaytimes = zeros(num_neurons,1);
nic_AoCs = zeros(num_neurons,1);
nic_startcalc = zeros(num_neurons,2);
nic_stopcalc = zeros(num_neurons,2);

kcl_max_dffs = zeros(num_neurons,1);
kcl_decaytimes = zeros(num_neurons,1);
kcl_AoCs = zeros(num_neurons,1);
kcl_startcalc = zeros(num_neurons,2);
kcl_stopcalc = zeros(num_neurons,2);


%defining text strings for row/column labels in output files
%--- testing single output excel file------
usublabel = cellstr('Substance');
uvertgap = cellstr('--------------------');
uneuron = cellstr('Neuron');
umaxpkstr = cellstr('Peak dF/F');
udecstr = cellstr('Decay times (s)');
uAoCstr = cellstr('Area under Curve (dF/F)');
ustartstr1 = cellstr('Time at Peak (s)');
ustartstr2 = cellstr('dF/F at Peak');
startnames = [ustartstr1 ustartstr2];
ustopstr1 = cellstr('Time to Settle (s)');
ustopstr2 = cellstr('dF/F at Settle (s)');
stopnames = [ustopstr1 ustopstr2];

%creating time vector from fps and length of files
t_length = (1/fps)*length(sg_combined.F_raw(1,:));
t = (0:(1/fps):t_length-(1/fps))';


%looping through each neuron that is detected and output by detect software
for i = 1:1:num_neurons
    cancel = 0;
    %individual dF/F load and normalizing shift
    rsg21raw = sg_combined.F_raw(i,:)';
    rsg21 = rsg21raw - min(rsg21raw);

    %lowpass filter of signal, 4=4th order, fc set in main code 
    [b,a] = butter(4,fc,'low');
    
    %filter coefficients b and a are fed into filter command to clean
    %signal
    filt_rsg21_1raw = filtfilt(b,a,rsg21);
    filt_rsg21_1 = filt_rsg21_1raw - min(filt_rsg21_1raw);

    %identifies all peaks and respective locations until kcl injection
    [rsgabspks,rsgabslocs_adj] = findpeaks(filt_rsg21_1(nic_bound:kcl_bound));
    rsgabslocs = int32(rsgabslocs_adj)+nic_bound;

    %clears any peaks that are less than 90% of the global peak of data
    rsgpks = rsgabspks(rsgabspks > max_tolerance*max(filt_rsg21_1(nic_bound:kcl_bound)));

    %clears the locations from same data set to equate size of vectors
    rsglocs = rsgabslocs(rsgabspks > max_tolerance*max(filt_rsg21_1(nic_bound:kcl_bound)));

    %concatenates the arrays, focused on first point where calculation
    %starts
    rsgmaxes = [t(rsglocs) rsgpks];
    
    %defines baseline of data
    baseline = mean(filt_rsg21_1(1:nic_bound));
    
    %if rsglocs is empty or the max peak is less than 2 x baseline, dataset
    %is classified as no meaningful peaks and does not perform calculations
    if isempty(rsglocs) || max(rsgpks) < insig_thresh 
        nic_max_dffs(i,1) = 0;
        nic_decaytimes(i,1) = 0;
        nic_AoCs(i,1) = 0;
        nic_startcalc(i,1) = 0;
        nic_startcalc(i,2) = 0;
        nic_stopcalc(i,1) = 0;
        nic_stopcalc(i,2) = 0;
        
    %when significant peak is detected...    
    else
        %sets end of calc threshold to 85% decay from max value
        rsgthresh = (rsgmaxes(1,2)-baseline)*low_scale+baseline;

        %defines time and amplitude array from detected max to end of file
        t_seg = t(rsglocs(1):end);
        filtseg = filt_rsg21_1(rsglocs(1):end); 

        %creates arrary that is cut off when amplitude goes below rsgthresh
        t_below = t_seg(filtseg < rsgthresh);
        
        %meausure until
        if isempty(t_below)
            cancel = 1;
            decay85 = 0; 
            aoc = 0;
            %fprintf('Neuron %i does not fall below 85%% decay threshold, no calculation performed.\n', i)
        else
            %takes array from max to 85% decay and performs subtraction for
            %decay time and integration for AuC; 
            decay85 = t_below(1) - rsgmaxes(1,1); 
            %converts time back to indexing value from max value until first
            %value when signal crosses rsgthresh; subtract baseline to avoid
            %inflation of values 
            aoc = trapz((filt_rsg21_1(int32(rsgmaxes(1,1)*fps):int32(t_below(1)*fps))-baseline));
        end
            
        %moving calucated values into proper arrays
        if cancel == 0
            nic_max_dffs(i,1) = max(filt_rsg21_1(nic_bound:kcl_bound));
        else
            nic_max_dffs(i,1) = 0;
        end
        nic_decaytimes(i,1) = decay85;
        nic_AoCs(i,1) = aoc;
        nic_startcalc(i,1) = rsgmaxes(1,1);
        nic_startcalc(i,2) = rsgmaxes(1,2);
        if ~isempty(t_below)
            nic_stopcalc(i,1) = t_below(1);
            nic_stopcalc(i,2) = (rsgmaxes(1,2)-baseline)*low_scale+baseline;
        else
            nic_stopcalc(i,1) = 0;
            nic_stopcalc(i,2) = 0;
        end
        
    end



    




%------KCl calculation here-------------------
    %identifies all peaks and respective locations until kcl injection
    [kcl_rsgabspks,kclrsgabslocs] = findpeaks(filt_rsg21_1(kcl_bound:end));
    kcl_rsgabslocs = kclrsgabslocs+double(kcl_bound);

    %clears any peaks that are less than 90% of the global peak of data
    kcl_rsgpks = kcl_rsgabspks(kcl_rsgabspks > max_tolerance*max(filt_rsg21_1(kcl_bound:end)));

    %clears the locations from same data set to equate size of vectors
    kcl_rsglocs = kcl_rsgabslocs(kcl_rsgabspks > max_tolerance*max(filt_rsg21_1(kcl_bound:end)));

    %concatenates the arrays, focused on first point where calculation
    %starts
    kcl_rsgmaxes = [t(kcl_rsglocs) kcl_rsgpks];
    
    kcl_baseline = mean(filt_rsg21_1(kcl_bound:kcl_bound+10));
    %if rsglocs is empty or the max peak is less tahn 2 x baseline, dataset
    %is classified as no meaningful peaks and does not perform calculations
    if isempty(kcl_rsglocs)
        kcl_max_dffs(i,1) = 0;
        kcl_decaytimes(i,1) = 0;
        kcl_AoCs(i,1) = 0;
        kcl_startcalc(i,1) = 0;
        kcl_startcalc(i,2) = 0;
        kcl_stopcalc(i,1) = 0;
        kcl_stopcalc(i,2) = 0;
        
    %when significant peak is detected...    
    else
        %sets end of calc threshold to 85% decay from max value
        kcl_rsgthresh = (kcl_rsgmaxes(1,2)-kcl_baseline)*low_scale+kcl_baseline;

        %defines time and amplitude array from detected max to end of file
        kcl_t_seg = t(kcl_rsglocs(1):end);
        kcl_filtseg = filt_rsg21_1(kcl_rsglocs(1):end); 

        %creates arrary that is cut off when amplitude goes below rsgthresh
        kcl_t_below = kcl_t_seg(kcl_filtseg < kcl_rsgthresh);
        
        %meausure until
        if isempty(kcl_t_below)
            decay85 = t(end) - kcl_rsgmaxes(1,1); 
            aoc = trapz((filt_rsg21_1(int32(kcl_rsgmaxes(1,1)*fps):int32(t(end)*fps))));
        else
            %takes array from max to 85% decay and performs subtraction for
            %decay time and integration for AuC; 
            decay85 = kcl_t_below(1) - kcl_rsgmaxes(1,1); 
            %converts time back to indexing value from max value until first
            %value when signal crosses rsgthresh; subtract baseline to avoid
            %inflation of values 
            aoc = trapz((filt_rsg21_1(int32(kcl_rsgmaxes(1,1)*fps):int32(kcl_t_below(1)*fps))));
        end
            
        %moving calucated values into proper arrays
        kcl_max_dffs(i,1) = max(filt_rsg21_1(kcl_bound:end));
        kcl_decaytimes(i,1) = decay85;
        kcl_AoCs(i,1) = aoc;
        kcl_startcalc(i,1) = kcl_rsgmaxes(1,1);
        kcl_startcalc(i,2) = kcl_rsgmaxes(1,2);
        if ~isempty(kcl_t_below)
            kcl_stopcalc(i,1) = kcl_t_below(1);
            kcl_stopcalc(i,2) = (kcl_rsgmaxes(1,2))*low_scale;
        else
            kcl_stopcalc(i,1) = t(end);
            kcl_stopcalc(i,2) = filt_rsg21_1(end);
        end

    end
%-----------------------------------------------------------

    %plotting filtered wave and calculation lines
          
    
    if isempty(rsglocs)  || max(rsgpks) < insig_thresh || isempty(t_below)
        fprintf('Neuron %i not plotted due to incomplete response.\n', i)
        removed(i) = i;
        continue 
            
    %if significant peaks were registered, calculation lines are plotted on
    %graph and proper titles are provided 
    else
        figure     
        hold on 
        plot(t,rsg21- min(filt_rsg21_1raw),'c')
        plot(t,filt_rsg21_1,'b',LineWidth=3) 
       
        xl = [0, kcl_bound/fps];
        yl = [rsgmaxes(1,2), rsgmaxes(1,2)];
        plot(xl,yl,'g--')
        lowbound = (rsgmaxes(1,2)-baseline)*low_scale+baseline;
        yl2 = [lowbound,lowbound];
        plot(xl,yl2,'r--')
        xline(rsgmaxes(1,1),'g--')
        xline(t_below(1),'r--')
        xline(kcl_bound/fps,'k', LineWidth=5)
        %test lines for spontaneous activity reading
        
        %------------
        if ~isempty(kcl_rsgmaxes)
            kclxl = [kcl_bound/fps, t(end)];
            kclyl = [kcl_rsgmaxes(1,2), kcl_rsgmaxes(1,2)];
            plot(kclxl,kclyl,'g--')
            kcl_lowbound = (kcl_rsgmaxes(1,2))*low_scale;
            kclyl2 = [kcl_lowbound,kcl_lowbound];
            plot(kclxl,kclyl2,'r--')
            xline(kcl_rsgmaxes(1,1),'g--')
            if ~isempty(kcl_t_below)
                xline(kcl_t_below(1),'r--')
            end    
        end
        hold off    
        titstrneuron = sprintf(' Neuron %0.f',i);
        titstr = [filenamestr, titstrneuron];
        title(titstr)
        xlabel('Time (s)')
        ylabel('dF/F')
        xlim([0 t(end)])
        ylim([0 max(maxes)+0.5])       
        legend('Raw dF/F','Filtered dF/F')
        saveas(gcf,[folderpath '\' titstr '.png']);
    end
    
end

%creating legend for combined figure and looping back thorugh code to plot;
%same plotting pipeline as original loop
legendstr = strings(size(num_neurons));
figure
hold on
for i = 1:1:num_neurons
    if removed(i) == 0
        rsg21raw = sg_combined.F_raw(i,:)';
        rsg21 = rsg21raw - min(rsg21raw);
        [b,a] = butter(4,fc,'low');
        filt_rsg21_1raw = filtfilt(b,a,rsg21);
        filt_rsg21_1 = filt_rsg21_1raw - min(filt_rsg21_1raw);
        plot(t,filt_rsg21_1,LineWidth=2)
        xlim([0 t(end)])
        
        ylim([0 max(maxes)+0.5])
        xlabel('Time (s)')
        ylabel('dF/F')
        
        legendstr(i) = sprintf('Neuron %.0f',i);
    else
        legendstr(i) = sprintf('Neuron %.0f not plotted',i);
        rsg21raw = sg_combined.F_raw(i,:)';
        rsg21 = rsg21raw - min(rsg21raw);
        plot(t,rsg21,'Color',[1 1 1 0])
    end
end
comstr = [filenamestr ': All Neurons'];
title(comstr)
xline(kcl_bound/fps ,'k', LineWidth=5)
hold off
legend(legendstr,'Location','bestoutside');

saveas(gcf,[folderpath '\' filenamestr ' All Neurons' '.png']);


for k = 1:1:num_neurons
    if removed(k) ~= 0
        filtset(k,:) = NaN;
    end
end

meansig = mean(filtset,'omitnan');
stdsig = std(filtset,'omitnan');
upsig = meansig+stdsig;
downsig = meansig-stdsig;


if num_neurons > 1
    figure
    hold on
    for x=1:num_neurons
        if removed(x) == 0
            plot(t,filtset(x,:),'Color',[1 0 0 0.5])
        end
    end
    shade(t,upsig,t,downsig,'FillType',[1 2;2 1],'Color',[0 0.25 1 0.8])
    plot(t,meansig,'b',LineWidth=2)
    hold off
    statsstr = [filenamestr ': Mean +/- 1 STD Overlay'];
    title(statsstr)
    xlabel('Time (s)')
    ylabel('dF/F')
    xlim([0 t(end)])   
    ylim([0 max(maxes)+0.5])
    saveas(gcf,[folderpath '\' filenamestr ' Mean+StDev Figure' '.png']);
end

%output arrays named and defined for main code to output to files

nicsublabelname = cell2table(nicnamecol,'VariableNames',usublabel);
kclsublabelname = cell2table(kclnamecol,'VariableNames',usublabel);
vertgapname = array2table(vertgap, 'VariableNames', uvertgap);
neuronlist = cell2table(neuroncol, 'VariableNames',uneuron);
nic_maxdff = array2table(nic_max_dffs, 'VariableNames',umaxpkstr);
nic_decaytime = array2table(nic_decaytimes,'VariableNames',udecstr);
nic_AoC = array2table(nic_AoCs,'VariableNames',uAoCstr);
nic_strtcalc = array2table(nic_startcalc,'VariableNames',startnames);
nic_stpcalc = array2table(nic_stopcalc,'VariableNames',stopnames);

kcl_maxdff = array2table(kcl_max_dffs, 'VariableNames',umaxpkstr);
kcl_decaytime = array2table(kcl_decaytimes,'VariableNames',udecstr);
kcl_AoC = array2table(kcl_AoCs,'VariableNames',uAoCstr);
kcl_strtcalc = array2table(kcl_startcalc,'VariableNames',startnames);
kcl_stpcalc = array2table(kcl_stopcalc,'VariableNames',stopnames);

nic_table = horzcat(nicsublabelname, neuronlist, nic_maxdff, nic_decaytime, nic_AoC, vertgapname, nic_strtcalc, nic_stpcalc);
kcl_table = horzcat(kclsublabelname, neuronlist, kcl_maxdff, kcl_decaytime, kcl_AoC, vertgapname, kcl_strtcalc, kcl_stpcalc);
gap = table('Size',[3,size(nic_table,2)], 'VariableTypes', nic_table.Properties.VariableTypes, 'VariableNames',nic_table.Properties.VariableNames);

trace_table = vertcat(nic_table, gap, kcl_table);
end

