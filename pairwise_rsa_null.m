function [results_matrix] = pairwise_rsa_null(test_matrix, training_matrix, niters)

%% Prep some basic values

% Number of stimulus classes to compare
number_classes = size(test_matrix,1);

% If the number of iterations to test is not specified, it will be all
% possible permutations of the labels up to 10 conditions. Otherwise, it
% will run 10! permutations of the labels.

if ~exist('niters','var'),
    niters = factorial(min([num_classes,10]));
end

% Save niters-many possible permutations in a matrix called samples
if number_classes<=10,
    samples = perms(1:num_classes);
    samples = samples(randsample(size(samples,1),niters),:);
else
    samples = cell2mat(arrayfun(@(x)randperm(number_classes),(1:niters)','UniformOutput',0));
end


% Number of comparisons to make (same as number of unique off-diagnoal
% values in the test or training matrix)
number_comparisons = (numel(test_matrix)-number_classes)/2;

% Generate a matrix to contain the resulst of all the comparisons for each
% iteration of the permutation search
results_matrix = nan(niters,number_comparisons);

if sum(isnan(test_matrix(:)))<numel(test_matrix),
    
    %% Loop through the iterations and look at the results
    tic;
    for iter = 1:niters,
        
        % Generate a list of every pairwise comparison and the results array
        list_of_comparisons = combnk([1:number_classes],2);
        results_of_comparisons = nan(number_comparisons,1);
        
        % The model matrix is permuted according to this iteration's sample
        training_matrix_i = training_matrix(samples(iter,:),samples(iter,:));
        
        
        %% Loop through all comparisons and test
        for this_comp = 1:number_comparisons
            
            test_classes = list_of_comparisons(this_comp,:);
            
            % Define which classes to keep in the RSA structures (non-test classes)
            non_test_classes = logical(ones(number_classes,1));
            non_test_classes(test_classes) = 0;
            
            % Rebuild the RSA matrices without the test classes in them
            test_matrix_temp = test_matrix(non_test_classes,:);
            training_matrix_temp = training_matrix_i(non_test_classes,:);
            
            % Extract the two vectors of interest from test and train matrices
            testA = test_matrix_temp(:,test_classes(1));
            testB = test_matrix_temp(:,test_classes(2));
            trainA = training_matrix_temp(:,test_classes(1));
            trainB = training_matrix_temp(:,test_classes(2));
            
            % Compare the correlations
            if (corr(testA,trainA)+corr(testB,trainB) > corr(testB,trainA)+corr(testA,trainB)),
                results_of_comparisons(this_comp) = 1;
            else
                results_of_comparisons(this_comp) = 0;
            end
            
        end
        results_matrix(iter,:) = results_of_comparisons';
        
        if ~mod(iter,500),
            fprintf('Iteration: %s...', num2str(iter));
            toc;
        end
    end
    
end
