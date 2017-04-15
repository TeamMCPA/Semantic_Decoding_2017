function subj_data = extract_vecs_1subj(file_path, cond_name_list, scan_window, exclude_by)


if ~exist('exclude_by','var'),
    fprintf('WARNING: No artifact exclusion for these data.\n');
    exclude_by = 'none';
else
    if ~(strcmp(exclude_by,'channel') || strcmp(exclude_by,'trial') || strcmp(exclude_by,'obsv') || strcmp(exclude_by,'none')),
        fprintf('ARTIFACT EXCLUSION METHOD NOT RECOGNIZED. Defaulting to "obsv".\n');
        exclude_by = 'obsv';
    end
end

%% Load in Windowed Data
if ischar(file_path),
    file_path = cellstr(file_path);
end

arrays=cell(length(file_path),1);
for filenum = 1:length(file_path),
    % Load the NIRS data from Probe Set #filenum
    load(file_path{filenum},'-mat')
    arrays{filenum} = squeeze(procResult.dc(:,1,:));
end
full_array = horzcat(arrays{:});

if ~strcmp(exclude_by,'none'),
    num_sds = 5; % 5 sd's for artifact
    window_size = 10; % 1 sec window around artefact
    fprintf('Removing signal artifacts (>%g SDs) in a %0.1f sec window.\n',num_sds,window_size/10);
    full_array = scrub_artifacts(full_array,num_sds,window_size); 
end

% Save the event onsets from the aux timeseries data
% cond_name_list = {'kitty' 'bunny' 'dog' 'bear' 'foot' 'hand' 'mouth' 'nose'};
condition = struct('name',cond_name_list);
for cond_num = 1:length(cond_name_list),
    condition(cond_num).onsets = find(aux(:,cond_num+1));
    condition(cond_num).onsets = condition(cond_num).onsets(1:2:end);
end

% Store windowed timeseries data for each condition

scan_onset = scan_window(1);
scan_offset = scan_window(2);
window_length = scan_offset-scan_onset+1;

for cond_num = 1:length(cond_name_list),
    % Create a 3-d matrix ( TIME x CHAN x EVENT )
    win_data = nan(window_length,size(full_array,2),length(condition(cond_num).onsets));
    
    for i = 1:length(condition(cond_num).onsets),
        
        % Set the onset and offset for this trial (align scan window to the
        % trial onset marker)
        this_onset = condition(cond_num).onsets(i) + scan_onset;
        this_offset = condition(cond_num).onsets(i) + scan_offset;
        
        % In case the end of the data comes before the trial offset
        if this_offset > size(full_array,1),
            fprintf('Offset #%g at %g exceeds end of data (at %g).\n',i,this_offset,size(full_array,1));
            this_offset = min(this_offset,size(full_array,1));
            fprintf('Using curtailed window [%g %g]\n\n',this_onset,this_offset);
        end
        time_series_length = this_offset-this_onset+1;
        
        % Extract the windowed data of interest for this trial
        win_data(1:time_series_length,:,i) = full_array(this_onset:this_offset,:) ;
        
        % Re-zero the windowed data based on the onset value. Importantly,
        % there is a decision to make here about using the actual trial
        % onset or the onset of the scan window. Since the scan window is a
        % bit arbitrary, I'm concerned we might cut off an actual rise in
        % the response prior to onset. The principled approach, IMO is to
        % use the trial onset which might result in a set of data all above
        % zero, but that is representative of the real response.
        win_data(:,:,i) = win_data(:,:,i) - repmat(full_array(condition(cond_num).onsets(i),:),window_length,1);
        
%         % z score values based on 1 sec before onset
%         win_mean = nanmean(full_array(condition(cond_num).onsets(i)-10:condition(cond_num).onsets(i),:));
%         win_std = nanstd(full_array(condition(cond_num).onsets(i)-10:condition(cond_num).onsets(i),:));
%         win_data(:,:,i) = (win_data(:,:,i) - repmat(win_mean,window_length,1)) ./ repmat(win_std,window_length,1);
    end
    
    % Save the windowed data for all trials of this condition to the struct
    condition(cond_num).windowed_data = win_data;
    
    % Identify events that contain artifacts
    contains_artifact = squeeze(logical(sum(isnan(win_data),1)));
    active_chans = ones(size(win_data,2),1);
    active_chans(sum(contains_artifact,2)==size(contains_artifact,2)) = 0;
    good_trials = ~logical(sum(contains_artifact(logical(active_chans),:),1));
    condition(cond_num).good_trials = good_trials;

    
    switch exclude_by
        
        case 'trial'
            
            % This scrub discards the whole trial if any channel is bad.
            % Average across all trials to get the average time series for this
            % condition in each channel
            condition(cond_num).window_averages = nanmean(win_data(:,:,good_trials),3);
            
        case 'channel'
            % This scrub only ignores the bad channel(s) and keeps the trial
            % Average across all trials to get the average time series for this
            % condition in each channel
            for event = 1:size(win_data,3),
                bad_chans = logical(sum(isnan(win_data(:,:,event))));
                win_data(:,bad_chans,event) = NaN;
            end
            condition(cond_num).window_averages = nanmean(win_data,3);
            
        case 'obsv'
            % This scrub will delete the NaN values only, but still use the
            % remaining data in the trial (even for artifact channel)
            condition(cond_num).window_averages = nanmean(win_data,3);
            
        case 'none'
            % This scrub would delete NaN values if they occurred, but no
            % NaN-masking of artifacts has been performed in this case.
            condition(cond_num).window_averages = nanmean(win_data,3);
            
            
    end
    
end

%% Compute some sort of summary statistic for the multi-channel patterns

for cond_num = 1:length(cond_name_list), 
    
    % Better to exclude the whole trial when computing ts_average than to
    % try to just exclude individual channels with artifacts in them
    condition(cond_num).ts_average = nanmean(condition(cond_num).window_averages)';
    %condition(cond_num).ts_average = nanmean(squeeze(mean(condition(cond_num).windowed_data)),2);
    
    % Don't know about second derivatives, but sure, maybe this will work.
    condition(cond_num).second_deriv = nanmean(diff(diff(condition(cond_num).window_averages)))';
end

subj_data = condition;

%save(['subjdata/' Probe1_path(7:(end-16)) '.mat'])

