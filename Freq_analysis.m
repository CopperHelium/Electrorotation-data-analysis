rate = 10000; 
F_state = UntitledVoltage_0.Data; % Shows whether the field was ON or OFF
X_signal = UntitledVoltage_1.Data; % Contains X-Photomultiplier signal
Y_signal = UntitledVoltage_2.Data; % Contains Y-Photomultiplier signal

%% 

Ob = fieldnames(ci); 
Props = fieldnames(ci.(Ob{end}));

i = 1;
while i <= size(Props,1)
    if isstruct(ci.(Ob{end}).(Props{i}))
        if strcmp(ci.(Ob{end}).(Props{i}).name,'wf_samples') 
            % number of data recorded per second
            rate = ci.(Ob{end}).(Props{i}).value; 
            break
        end
    end
    i = i+1;
end
