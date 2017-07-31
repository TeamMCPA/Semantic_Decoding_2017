function [results_of_comparisons, list_of_comparisons] = pairwise_rsa_test(test_matrix, training_matrix)
% PAIRWISE_RSA_TEST Perform all pairwise (two-alternative forced choice)
% comparisons between two similarity structures.
%   ACC = PAIRWISE_RSA_TEST( M1, M2 ) returns the results (binary) of each
%   pairwise comparison in a vector.
%
%   M1 and M2 must be square, symmetric, and the same size, but they may be
%   populated with any distance values (Pearson R, Fisher R-to-Z, euclidean
%   distance, etc).
%
%   For n-by-n similarity matrices, all pairs of rows (all n-pick-2) are
%   tested using: combnk(1:n,2). List of comparisons (COMP) also available:
%   [ ACC, COMP ] = PAIRWISE_RSA_TEST( M1, M2 )

%% Sanity Check
if sum(isnan(test_matrix(:))) || sum(isnan(training_matrix(:)))
    disp('One or both input matrices contains all NaN values. I quit!');
    return
end


%% Prep some basic values

% Number of stimulus classes to compare, based on the number of rows
number_classes = size(test_matrix,1);

% Generate a list of every pairwise comparison and the results array
list_of_comparisons = combnk([1:number_classes],2);
number_of_comparisons = size(list_of_comparisons,1);
results_of_comparisons = nan(number_of_comparisons,1);

%% Loop through all comparisons and test
for this_comp = 1:number_of_comparisons
    
    test_classes = list_of_comparisons(this_comp,:);
    
    % Define which classes to keep in the RSA structures (non-test classes)
    non_test_classes = logical(ones(number_classes,1));
    non_test_classes(test_classes) = 0;
    
    % Rebuild the RSA matrices with the test classes' rows deleted
    test_matrix_temp = test_matrix(non_test_classes,:);
    training_matrix_temp = training_matrix(non_test_classes,:);
    
    % Extract the two column vectors for the test classes from each matrix
    testA = test_matrix_temp(:,test_classes(1));
    testB = test_matrix_temp(:,test_classes(2));
    trainA = training_matrix_temp(:,test_classes(1));
    trainB = training_matrix_temp(:,test_classes(2));
    
    % Compare the correlations for A-A,B-B vs. A-B,B-A
    correct_labels = corr(testA,trainA,'rows','pairwise') + corr(testB,trainB,'rows','pairwise');
    incorrect_labels = corr(testB,trainA,'rows','pairwise') + corr(testA,trainB,'rows','pairwise');
    
    % Test the accuracy, whether correct label correlations were greater
    % than the incorrect label correlations
    results_of_comparisons(this_comp) = correct_labels > incorrect_labels;
    
end