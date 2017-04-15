function accuracies = pairwise_rsa_leaveoneout_null(test_subjs_mat, niters)

subjs = 1:size(test_subjs_mat,3);
nsubjs = length(subjs);
number_classes = size(test_subjs_mat,1);

% If the number of iterations to test is not specified, it will be all
% possible permutations of the matrix.
if ~exist('niters','var'),
    % Save all possible permutations in a matrix called samples
    samples = perms(1:number_classes);
    niters = size(samples,1);
else
    % Save niters-many possible permutations in a matrix called samples
    samples = perms(1:number_classes);
    samples = samples(randsample(size(samples,1),niters),:);
end

acc = nan(nsubjs,1);
accuracies = nan(niters,nsubjs);

tic;
for iter = 1:niters
    
    for this_sub = subjs,
        
        group_mat = nanmean(test_subjs_mat(:,:,subjs(subjs~=this_sub)),3);
        subj_mat = test_subjs_mat(:,:,this_sub);
        subj_mat = subj_mat(samples(iter,:),samples(iter,:));
        
        acc(this_sub) = mean(pairwise_rsa_test(subj_mat,group_mat));
        
    end
    
    accuracies(iter,:) = acc';
    
    if ~mod(iter,500),
        fprintf('Iteration: %s...', num2str(iter));
        toc;
    end
    
end