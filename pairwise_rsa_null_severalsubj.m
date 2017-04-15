function accuracies = pairwise_rsa_null_severalsubj(test_subjs_mat, model_mat, niters)

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
    
    model_mat_i = model_mat(samples(iter,:),samples(iter,:));

    for this_sub = subjs,
        
        subj_mat = test_subjs_mat(:,:,this_sub);
        
        acc(this_sub) = mean(pairwise_rsa_test(subj_mat,model_mat_i));
        
    end
    
    accuracies(iter,:) = acc';
    
    if ~mod(iter,500),
        fprintf('Iteration: %s...', num2str(iter));
        toc;
    end
end