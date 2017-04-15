function full_array_scrub = scrub_artifacts(full_array,num_std,window_size)

% copy the full array
% scrubbing will just fill in NaNs where the outliers were
full_array_scrub = full_array;

% Loop through and evaluate each channel separately
for chan = 1:size(full_array,2),
    
    % unscrubbed time series data
    signal = full_array(:,chan);

    % calculate some statistics
    signal_mean = nanmean(signal);
    signal_std = nanstd(signal);
    signal_thresh_hi = signal_mean + num_std * signal_std;
    signal_thresh_lo = signal_mean - num_std * signal_std;
    
    % outliers are any observation above or below the # of sd's
    outlier = (signal>signal_thresh_hi | signal<signal_thresh_lo);
    
    % change them to NaNs
    signal(outlier) = NaN;
    
    % save to the output array
    full_array_scrub(:,chan) = signal;
    
end

if window_size > 1,
    % Define a logical mask for the full_array_scrub that marks the locations
    % of artefacts. This mask will NaN-out windows around any artifacts
    % according to window_size.
    artefact_mask = isnan(filter2(ones(window_size,1),full_array_scrub));
    full_array_scrub(artefact_mask) = NaN;
    
    % If the window_size is zero, it will remove the original marked
    % artefacts, so this procedure should be skipped entirely unless
    % window_size is >1.
    
end