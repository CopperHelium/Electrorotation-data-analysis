clc;
close all
clear
warning('off', 'MATLAB:colon:nonIntegerIndex')

%Ask the user to select folder containing the .mat data files to be
% analyzed for freq. Use current working directory as a starting point.

start_path=pwd;
folder_name = uigetdir(start_path,'Select folder containing .mat files to be analyzed:');
cd(folder_name)
%----------------------------------------------------------------

% Collect a list of all the files in the current directory
Files=dir('*s.mat');

% Measure how long the list of files is
lengthFiles=length(Files);



for kk = 1:lengthFiles
    clc;
    clearvars -except kk Files lengthFiles folder_name
    
    file_name = Files(kk).name;
    path_name = strcat(folder_name,'\',file_name);
    
    %Load the mat file containing the data. Eventually this should be in a loop
    %going through all files in a cerain folder(s).
    load(path_name)
    
    X_signal = UntitledVoltage_0.Data; % Contains X-Photomultiplier signal
    Y_signal = UntitledVoltage_1.Data; % Contains Y-Photomultiplier signal
    fs = UntitledVoltage_2.Property.wf_samples; % number of samples per sec
    F_state = UntitledVoltage_2.Data; % Shows whether the field was ON or OFF
    F_state = (F_state > 2); % convert to binary (0/1)
    % Time = (10^-4)*(1:size(X_signal))';
    
    Switch = diff(F_state); % find if the field was switched
    OnOff = find(Switch); % Indices of the switching data points
    
    %Define pairs in sample numbers describing the times when the field was on
    %and off. The first column is when the event startes and the second when it
    %ended. Frequency calculated between these two.
    
    ON_pairs = [OnOff(1:2:end) OnOff(2:2:end)];
    OFF_pairs = [OnOff(2:2:end-2) OnOff(3:2:end)];
    
    ON_cnt = (ON_pairs(3,2)-ON_pairs(3,1)); % # samples in each ON event
    OFF_cnt = (OFF_pairs(3,2)-OFF_pairs(3,1)); % # samples in each OFF event
    ON_dur = ON_cnt/fs; % actual duration of each ON event
    OFF_dur = OFF_cnt/fs; % actual duration of each OFF event
    OnOff_duration = OFF_dur + ON_dur; % Total duration of each cycle
    
    % The matrix above missed the last OFF event immediately following the last
    % ON event and before the recovery starts (two seconds after the last ON).
    % Adding that bit to the Off_pairs here.
    rec_flag = 0;
    if size(F_state,1)-ON_pairs(end,2)>OFF_cnt && F_state(end)==0
        OFF_pairs = vertcat([1 OnOff(1)],OFF_pairs,[OnOff(end) OnOff(end)+(OFF_pairs(3,2)-OFF_pairs(3,1))]);
        rec_flag = 1;
    elseif size(F_state,1)-ON_pairs(end,2)<OFF_cnt && F_state(end)==0
        OFF_pairs = vertcat([1 OnOff(1)],OFF_pairs);
    else
        OFF_pairs = vertcat([1 OnOff(1)],OFF_pairs);
    end
    
    %Frequency calculation
    
    FRE_ON = zeros(size(ON_pairs,1),1);
    FRE_OFF = zeros(size(OFF_pairs,1),1);
    
    T_ON = zeros(size(ON_pairs,1),1);
    T_OFF = zeros(size(OFF_pairs,1),1);
    
    %ON
    for ii = 1:size(ON_pairs,1)
        
        % apply a highpass filter on the ON data to prevent spurious low freqs from
        % affecting the calculation. The ON freq in proper experiments done above
        % ~15 degree C is expected to be well above ~100 Hz. Keep the cutoff fc
        % around 20-50
        
        fc = 40;
        [b,a] = butter(6,fc/(fs/2),'high');
        
        %Get the "chunk" of data for this computation
        X_chunk = X_signal(ON_pairs(ii,1):ON_pairs(ii,2));
        Y_chunk = Y_signal(ON_pairs(ii,1):ON_pairs(ii,2));
        
        
        X_chunk = filter(b,a,X_chunk); %apply the highpass filter
        Y_chunk = filter(b,a,Y_chunk); %apply the highpass filter
        
        
        fftx = fft(X_chunk - mean(X_chunk));
        ffty = fft(Y_chunk - mean(Y_chunk));
        
        P2x = abs(fftx/length(X_chunk));
        P1x = P2x(1:length(X_chunk)/2+1);
        P1x(2:end-1) = 2*P1x(2:end-1);
        
        P2y = abs(ffty/length(Y_chunk));
        P1y = P2y(1:length(Y_chunk)/2+1);
        P1y(2:end-1) = 2*P1y(2:end-1);
        
        f = fs*(0:(length(Y_chunk)/2))/length(Y_chunk);
        
        if ii==1
            [~,loc] = max(P1x+P1y);
        end
        
        [peakLoc, peakMag] = peakfinder(P1x+P1y);
        [Mag_d,I] = sort(peakMag,'descend');
        if length(peakLoc)>5
            Loc_d = peakLoc(I(1:5));
        else
            Loc_d = peakLoc(I);
        end
        loc_diff = abs(Loc_d - loc);
        [~,min_loc]= min(loc_diff);
        loc = Loc_d(min_loc);
        FRE_ON(ii) = f(loc);
        T_ON(ii) = (mean(ON_pairs(ii,:))-ON_pairs(1,1))/fs;
        % The other method to work almost equally well was pmtm (below).
        % [pxx2,f2] = pmtm(X_chunk - mean(X_chunk),1.25,[],rate);
        % [~,loc] = max(pxx2);
        % FRE_ON_pmtm = f2(loc)
        
        % If needed, the power spectrum can be checked using the code below
        % figure;
        % plot(f,pxx)
        % xlim([0 500])
        
    end
    
    %OFF
    
    for ii = 1:size(OFF_pairs,1) %ON_pairs and OFF_pairs are NOT the same size
        
        % apply a low filter on the OFF data to prevent spurious high freqs from
        % affecting the calculation. The OFF freq in proper experiments done above
        % ~15 degree C is expected to be well below ~50 Hz. Keep the cutoff fc
        % around 50-100 Hz
        
        fc = 50;
        [b,a] = butter(6,fc/(fs/2),'low');
        
        X_chunk = X_signal(OFF_pairs(ii,1):OFF_pairs(ii,2));
        Y_chunk = Y_signal(OFF_pairs(ii,1):OFF_pairs(ii,2));
        
        X_chunk = filter(b,a,X_chunk); %apply the lowpass filter
        Y_chunk = filter(b,a,Y_chunk); %apply the lowpass filter
        
        
        fftx = fft(X_chunk - mean(X_chunk));
        ffty = fft(Y_chunk - mean(Y_chunk));
        
        P2x = abs(fftx/length(X_chunk));
        P1x = P2x(1:length(X_chunk)/2+1);
        P1x(2:end-1) = 2*P1x(2:end-1);
        
        P2y = abs(ffty/length(Y_chunk));
        P1y = P2y(1:length(Y_chunk)/2+1);
        P1y(2:end-1) = 2*P1y(2:end-1);
        
        f = fs*(0:(length(Y_chunk)/2))/length(Y_chunk);
        
        if ii==1
            [~,loc] = max(P1x+P1y);
        end
        [peakLoc, peakMag] = peakfinder(P1x+P1y);
        [Mag_d,I] = sort(peakMag,'descend');
        if length(peakLoc)>5
            Loc_d = peakLoc(I(1:5));
        else
            Loc_d = peakLoc(I);
        end
        loc_diff = abs(Loc_d - loc);
        [~,min_loc]= min(loc_diff);
        loc = Loc_d(min_loc);
        FRE_OFF(ii) = f(loc);
        T_OFF(ii) = (mean(OFF_pairs(ii,:))-ON_pairs(1,1))/fs;
        
        % [pxx2,f2] = pmtm(X_chunk - mean(X_chunk),1.25,[],rate);
        % [~,loc] = max(pxx2);
        % FRE_Off_pmtm = f2(loc)
        
        % figure;
        % subplot(3,1,1)
        % plot(f,pxx)
        % xlim([0 20])
        
    end
    T_OFF(1)=0; % The first one is at t = 0 in the experiment.
    
    
    % REC
    if rec_flag == 1
        
        %Recovery starts when the final OFF event ends, until the end of sampling.
        %Divided that section of the data into chunks of the same size as OFF
        %events for frequency calculation. Then the same parameters as OFF can be
        %used for REC.
        REC_array = OFF_pairs(end,2):1*OFF_cnt:size(X_signal,1);
        REC_pairs = [REC_array(1:end-1)' REC_array(2:end)'];
        
        FRE_REC = zeros(size(REC_pairs,1),1);
        T_REC = zeros(size(REC_pairs,1),1);
        loc = find(f==FRE_OFF(end));
        for jj = 1:size(REC_pairs,1) %REC_pairs are of a different size
            
            % apply a low filter on the REC data to prevent spurious high freqs from
            % affecting the calculation. The REC freq in proper experiments done above
            % ~15 degree C is expected to be well below ~50 Hz. Keep the cutoff fc
            % around 50-100 Hz
            
            fc = 2*max(FRE_OFF);
            [b,a] = butter(6,fc/(fs/2),'low');
            
            X_chunk = X_signal(REC_pairs(jj,1):REC_pairs(jj,2));
            Y_chunk = Y_signal(REC_pairs(jj,1):REC_pairs(jj,2));
            
            X_chunk = filter(b,a,X_chunk); %apply the lowpass filter
            Y_chunk = filter(b,a,Y_chunk); %apply the lowpass filter
            
            
            [pxx,f] = periodogram(X_chunk - mean(X_chunk), [],[],fs);
            pyy = periodogram(Y_chunk - mean(Y_chunk), [],[],fs);
            if jj==1
                [~,loc] = max(pxx+pyy);
            end
            [peakLoc, peakMag] = peakfinder(pxx+pyy);
            [Mag_d,I] = sort(peakMag,'descend');
            if length(peakLoc)>5
                Loc_d = peakLoc(I(1:5));
            else
                Loc_d = peakLoc(I);
            end
            loc_diff = abs(Loc_d - loc);
            [~,min_loc]= min(loc_diff);
            loc = Loc_d(min_loc);
            FRE_REC(jj) = f(loc);
            T_REC(jj) = (mean(REC_pairs(jj,:))-ON_pairs(1,1))/fs;
            
            %         plot(f,pxx+pyy)
            %         xlim([0 50])
            %
            %         [pxx2,f2] = pmtm(X_chunk - mean(X_chunk),1.25,[],rate);
            %         [~,loc] = max(pxx2);
            %         FRE_Off_pmtm = f2(loc)
            %
            %         figure;
            %         plot(f,pxx)
            
            
        end
    end
    
    FRE_ON_mean = mean(FRE_ON);
    FRE_ON_std = std(FRE_ON);
    
    file_out = strcat(file_name(1:end-4),'_OUT.mat');
    if rec_flag == 1
%         save(file_out,'FRE_OFF','FRE_ON','FRE_REC','T_ON','T_OFF','T_REC','rec_flag')
    else
%         save(file_out,'FRE_OFF','FRE_ON','T_ON','T_OFF','rec_flag')
    end
    % Plot the ON frequency with time
    % figure;
    % plot(T_ON,FRE_ON,'--o')
    % ylim([0 max(FRE_ON)+ (max(FRE_ON) - min(FRE_ON))])
    % xlabel('Time (s)')
    % ylabel('Frequency (Hz)')
    
    % Plot the OFF frequency with time.
    
    figure;
    A1 = axes('position',[0.08 0.2 0.4 0.6]);
    A2 = axes('position',[0.56 0.2 0.4 0.6]);
    set(gcf,'CurrentAxes',A1)
    set(gca, 'FontSize', 20)
    plot(T_OFF,FRE_OFF,'-s','color','k')
    xlabel('Time (s)','FontSize', 20)
    ylabel('Frequency (Hz)','FontSize', 20)
    title('Dissociation')
    if rec_flag
        top = 2*ceil(0.5*max([FRE_OFF;FRE_REC]));
    else
        top = 2*ceil(0.5*max([FRE_OFF]));
    end
    ylim([0 top])
    set(gca, 'FontSize', 20)
    
    
    set(gcf,'CurrentAxes',A2)
    set(gca, 'FontSize', 20)
    if rec_flag
        plot(T_REC,FRE_REC,'-s','color','k')
    end
    xlabel('Time (s)','FontSize', 20)
    ylabel('Frequency (Hz)','FontSize', 20)
    title('Recovery')
    ylim([0 top])
    
    set(gcf,'PaperType','A4')
    set(gcf,'PaperUnits','centimeters')
    xSize = 28; ySize = 15;
    xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
    set(gcf,'PaperPosition',[xLeft yTop xSize ySize])
    set(gcf,'Position',[100 100 xSize*50 ySize*50])
    set(gcf,'PaperOrientation','landscape');
    
    axes('position',[0 0 1 1],'visible','off')
    text(0.17,0.75,strcat('Electrorotation speed = ',sprintf(' %0.0f',FRE_ON_mean),char(177),sprintf('%0.0f',FRE_ON_std)),'FontSize', 18)
    
    if rec_flag == 0
        text(0.69,0.50,'No recovery data','FontSize', 18)
    end
    axis([0 1 0 1])
    %     set(h, 'Position', [50, 200, 1400, 500]);
    plot_name = strcat(file_name(1:end-4),'.pdf');
%     saveas(gcf,plot_name)
end