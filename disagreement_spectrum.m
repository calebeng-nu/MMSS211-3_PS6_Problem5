%% Problem 5: Disagreement Spectrum in the Class Data
% This script solves Problem 5.
%
% It can be run all at once or section-by-section using %%.
% Run the sections in order if running section-by-section.
%
% Outputs:
%   1. A comprehensive report in the Command Window.
%   2. A saved text report called problem5_report.txt.
%   3. An optional histogram comparing random cuts to the sign split.

clear; clc;

report_filename = 'problem5_report.txt';

% Start fresh report
fid = fopen(report_filename, 'w');
if fid == -1
    error('Could not open problem5_report.txt for writing.');
end

fprintf(fid, '============================================================\n');
fprintf(fid, 'Problem 5: Disagreement Spectrum in the Class Data\n');
fprintf(fid, '============================================================\n\n');
fclose(fid);

fprintf('============================================================\n');
fprintf('Problem 5: Disagreement Spectrum in the Class Data\n');
fprintf('============================================================\n\n');

url = 'https://bengolub-economics.github.io/undergrad_networks_2026/shareable/info_professional.csv';
filename = 'info_professional.csv';

fprintf('Data URL:\n%s\n\n', url);

fid = fopen(report_filename, 'a');
fprintf(fid, 'Data URL:\n%s\n\n', url);
fclose(fid);

% Download file if needed
if ~isfile(filename)
    fprintf('File not found locally. Downloading %s...\n', filename);
    websave(filename, url);
    fprintf('Download complete.\n\n');

    fid = fopen(report_filename, 'a');
    fprintf(fid, 'File not found locally. Downloaded %s.\n\n', filename);
    fclose(fid);
else
    fprintf('Using local file: %s\n\n', filename);

    fid = fopen(report_filename, 'a');
    fprintf(fid, 'Using local file: %s\n\n', filename);
    fclose(fid);
end

% Read the matrix robustly.
% readmatrix is usually best for mostly numeric CSV files.
Raw = readmatrix(filename);

% Remove completely empty rows or columns if MATLAB imported headers as NaNs.
Raw = Raw(~all(isnan(Raw), 2), :);
Raw = Raw(:, ~all(isnan(Raw), 1));

% If there are still NaNs, try reading as a table and converting.
if any(isnan(Raw), 'all')
    T = readtable(filename, 'VariableNamingRule', 'preserve');
    Raw2 = table2array(T);

    if ~isnumeric(Raw2)
        Raw2 = str2double(string(Raw2));
    end

    Raw2 = Raw2(~all(isnan(Raw2), 2), :);
    Raw2 = Raw2(:, ~all(isnan(Raw2), 1));

    Raw = Raw2;
end

% If the matrix is not square, try common fixes:
% 1. Remove first column if it is an index column.
% 2. Remove first row if it is a header-like row.
M = Raw;

if size(M,1) ~= size(M,2)
    M_try = M(:, 2:end);
    if size(M_try,1) == size(M_try,2)
        M = M_try;
        fprintf('Detected and removed first column as an index column.\n\n');

        fid = fopen(report_filename, 'a');
        fprintf(fid, 'Detected and removed first column as an index column.\n\n');
        fclose(fid);
    end
end

if size(M,1) ~= size(M,2)
    M_try = M(2:end, :);
    if size(M_try,1) == size(M_try,2)
        M = M_try;
        fprintf('Detected and removed first row as a header-like row.\n\n');

        fid = fopen(report_filename, 'a');
        fprintf(fid, 'Detected and removed first row as a header-like row.\n\n');
        fclose(fid);
    end
end

if size(M,1) ~= size(M,2)
    error('Could not form a square matrix from the CSV. Matrix size is %d by %d.', size(M,1), size(M,2));
end

if any(isnan(M), 'all')
    error('Matrix still contains NaN values after cleanup. Check CSV formatting.');
end

n_original = size(M,1);

fprintf('Original weighted matrix size: %d by %d\n', size(M,1), size(M,2));
fprintf('Total positive entries in original matrix: %d\n\n', nnz(M > 0));

fid = fopen(report_filename, 'a');
fprintf(fid, 'Original weighted matrix size: %d by %d\n', size(M,1), size(M,2));
fprintf(fid, 'Total positive entries in original matrix: %d\n\n', nnz(M > 0));
fclose(fid);


%% Part 5(a): Strongly connected components and construction of W

fprintf('============================================================\n');
fprintf('Part 5(a): Strongly Connected Components and W\n');
fprintf('============================================================\n\n');

fid = fopen(report_filename, 'a');
fprintf(fid, '============================================================\n');
fprintf(fid, 'Part 5(a): Strongly Connected Components and W\n');
fprintf(fid, '============================================================\n\n');
fclose(fid);

% A_ij = 1 if the corresponding weighted entry is positive.
A = double(M > 0);

% Directed graph where edge i -> j exists if A(i,j) = 1.
G = digraph(A);

component_id = conncomp(G, 'Type', 'strong');
component_id = component_id(:);

num_components = max(component_id);
component_sizes = accumarray(component_id, 1);

[largest_size, largest_component] = max(component_sizes);

largest_indices = find(component_id == largest_component);
outside_indices = find(component_id ~= largest_component);

fprintf('Number of strongly connected components: %d\n', num_components);
fprintf('Size of largest strongly connected component: %d\n', largest_size);
fprintf('Original number of nodes: %d\n\n', n_original);

fprintf('All component sizes:\n');
disp(component_sizes');

fprintf('Indices in largest component:\n');
disp(largest_indices');

fprintf('Indices not in largest component:\n');
disp(outside_indices');

fid = fopen(report_filename, 'a');
fprintf(fid, 'Number of strongly connected components: %d\n', num_components);
fprintf(fid, 'Size of largest strongly connected component: %d\n', largest_size);
fprintf(fid, 'Original number of nodes: %d\n\n', n_original);

fprintf(fid, 'All component sizes:\n');
fprintf(fid, '%s\n\n', mat2str(component_sizes'));

fprintf(fid, 'Indices in largest component:\n');
fprintf(fid, '%s\n\n', mat2str(largest_indices'));

fprintf(fid, 'Indices not in largest component:\n');
fprintf(fid, '%s\n\n', mat2str(outside_indices'));
fclose(fid);

% Restrict original weighted matrix to largest SCC.
M_scc = M(largest_indices, largest_indices);

% Renormalize rows to sum to 1.
row_sums = sum(M_scc, 2);

if any(row_sums == 0)
    error('A row in the largest SCC has sum zero, so W cannot be row-normalized.');
end

W = M_scc ./ row_sums;

fprintf('\nConstructed DeGroot matrix W of size %d by %d.\n', size(W,1), size(W,2));
fprintf('Minimum row sum of W: %.12f\n', min(sum(W,2)));
fprintf('Maximum row sum of W: %.12f\n\n', max(sum(W,2)));

fid = fopen(report_filename, 'a');
fprintf(fid, 'Constructed DeGroot matrix W of size %d by %d.\n', size(W,1), size(W,2));
fprintf(fid, 'Minimum row sum of W: %.12f\n', min(sum(W,2)));
fprintf(fid, 'Maximum row sum of W: %.12f\n\n', max(sum(W,2)));
fclose(fid);


%% Part 5(b): Second eigenvalue and corresponding real right-hand eigenvector

fprintf('============================================================\n');
fprintf('Part 5(b): Second Eigenvalue and Right-Hand Eigenvector\n');
fprintf('============================================================\n\n');

fid = fopen(report_filename, 'a');
fprintf(fid, '============================================================\n');
fprintf(fid, 'Part 5(b): Second Eigenvalue and Right-Hand Eigenvector\n');
fprintf(fid, '============================================================\n\n');
fclose(fid);

% MATLAB's eig gives right eigenvectors:
% W * V_eig(:,j) = D_eig(j,j) * V_eig(:,j).
[V_eig, D_eig] = eig(W);
lambda = diag(D_eig);

% Sort eigenvalues by absolute value, largest first.
[~, order] = sort(abs(lambda), 'descend');

lambda_sorted = lambda(order);
V_sorted = V_eig(:, order);

lambda1 = lambda_sorted(1);
lambda2 = lambda_sorted(2);
sigma = V_sorted(:, 2);

% The prompt asks for a real right-hand eigenvector.
% In this data, the relevant eigenvector should be real up to numerical error.
if max(abs(imag(sigma))) < 1e-10
    sigma = real(sigma);
else
    warning('The selected second eigenvector has nontrivial imaginary parts. Taking real part for reporting.');
    sigma = real(sigma);
end

% Normalize to Euclidean norm 1.
sigma = sigma / norm(sigma);

% Sign convention: make the first nonzero entry positive.
% Eigenvectors are only defined up to multiplication by -1.
first_nonzero = find(abs(sigma) > 1e-10, 1);
if ~isempty(first_nonzero) && sigma(first_nonzero) < 0
    sigma = -sigma;
end

fprintf('Largest eigenvalue lambda_1:\n');
fprintf('  %.12f + %.12fi\n\n', real(lambda1), imag(lambda1));

fprintf('Second-largest-by-absolute-value eigenvalue lambda_2:\n');
fprintf('  %.12f + %.12fi\n\n', real(lambda2), imag(lambda2));

fprintf('Check Euclidean norm of sigma: %.12f\n\n', norm(sigma));

fid = fopen(report_filename, 'a');
fprintf(fid, 'Largest eigenvalue lambda_1:\n');
fprintf(fid, '  %.12f + %.12fi\n\n', real(lambda1), imag(lambda1));

fprintf(fid, 'Second-largest-by-absolute-value eigenvalue lambda_2:\n');
fprintf(fid, '  %.12f + %.12fi\n\n', real(lambda2), imag(lambda2));

fprintf(fid, 'Check Euclidean norm of sigma: %.12f\n\n', norm(sigma));
fclose(fid);

% Robust table construction.
largest_indices_col = largest_indices(:);
sigma_col = sigma(:);

if length(largest_indices_col) ~= length(sigma_col)
    error('Mismatch: largest_indices has length %d but sigma has length %d.', ...
        length(largest_indices_col), length(sigma_col));
end

sigma_table = table(largest_indices_col, sigma_col, ...
    'VariableNames', {'OriginalIndex', 'Sigma'});

fprintf('Sigma entries, using original student indices:\n');
disp(sigma_table);

fid = fopen(report_filename, 'a');
fprintf(fid, 'Sigma entries, using original student indices:\n');
for ii = 1:length(sigma_col)
    fprintf(fid, 'OriginalIndex %d: %.12f\n', largest_indices_col(ii), sigma_col(ii));
end
fprintf(fid, '\n');
fclose(fid);


%% Part 5(c): Divide nodes into left and right sets by sign of sigma

fprintf('============================================================\n');
fprintf('Part 5(c): Sign Split\n');
fprintf('============================================================\n\n');

fid = fopen(report_filename, 'a');
fprintf(fid, '============================================================\n');
fprintf(fid, 'Part 5(c): Sign Split\n');
fprintf(fid, '============================================================\n\n');
fclose(fid);

tol = 1e-10;

L_local = find(sigma < -tol);
R_local = find(sigma > tol);
Zero_local = find(abs(sigma) <= tol);

% Convert local SCC indices back to original indices.
L = largest_indices(L_local);
R = largest_indices(R_local);
Zero = largest_indices(Zero_local);

fprintf('Left set L = {i : sigma_i < 0}\n');
fprintf('Size of L: %d\n', numel(L));
disp(L');

fprintf('Right set R = {i : sigma_i > 0}\n');
fprintf('Size of R: %d\n', numel(R));
disp(R');

if ~isempty(Zero)
    fprintf('Nodes with sigma approximately zero:\n');
    disp(Zero');
else
    fprintf('No nodes have sigma approximately zero.\n');
end

fid = fopen(report_filename, 'a');
fprintf(fid, 'Left set L = {i : sigma_i < 0}\n');
fprintf(fid, 'Size of L: %d\n', numel(L));
fprintf(fid, '%s\n\n', mat2str(L'));

fprintf(fid, 'Right set R = {i : sigma_i > 0}\n');
fprintf(fid, 'Size of R: %d\n', numel(R));
fprintf(fid, '%s\n\n', mat2str(R'));

if ~isempty(Zero)
    fprintf(fid, 'Nodes with sigma approximately zero:\n');
    fprintf(fid, '%s\n\n', mat2str(Zero'));
else
    fprintf(fid, 'No nodes have sigma approximately zero.\n\n');
end
fclose(fid);


%% Part 5(d): Crossing weight for sign split

fprintf('============================================================\n');
fprintf('Part 5(d): Crossing Weight for the Sign Split\n');
fprintf('============================================================\n\n');

fid = fopen(report_filename, 'a');
fprintf(fid, '============================================================\n');
fprintf(fid, 'Part 5(d): Crossing Weight for the Sign Split\n');
fprintf(fid, '============================================================\n\n');
fclose(fid);

% C(L,R) = sum_{i in L, j in R} W_ij + sum_{i in R, j in L} W_ij.
C_sign = sum(W(L_local, R_local), 'all') + sum(W(R_local, L_local), 'all');

fprintf('Crossing weight for sign split C(L,R): %.12f\n\n', C_sign);

fid = fopen(report_filename, 'a');
fprintf(fid, 'Crossing weight for sign split C(L,R): %.12f\n\n', C_sign);
fclose(fid);


%% Part 5(d): Random partitions with same two set sizes

fprintf('============================================================\n');
fprintf('Part 5(d): Random Partition Comparison\n');
fprintf('============================================================\n\n');

fid = fopen(report_filename, 'a');
fprintf(fid, '============================================================\n');
fprintf(fid, 'Part 5(d): Random Partition Comparison\n');
fprintf(fid, '============================================================\n\n');
fclose(fid);

rng(1); % Reproducibility.

num_random = 10000; % Prompt asks for at least 1000.
m = size(W,1);

size_L = numel(L_local);
size_R = numel(R_local);

if size_L + size_R ~= m
    warning('Some sigma entries are approximately zero. Random partitions use only the sizes of L and R.');
end

random_crossings = zeros(num_random, 1);

for r = 1:num_random
    perm = randperm(m);

    L_rand = perm(1:size_L);
    R_rand = perm(size_L+1:size_L+size_R);

    random_crossings(r) = sum(W(L_rand, R_rand), 'all') + sum(W(R_rand, L_rand), 'all');
end

random_mean = mean(random_crossings);
random_sd = std(random_crossings);
random_min = min(random_crossings);
random_median = median(random_crossings);
random_p10 = prctile(random_crossings, 10);
random_p5 = prctile(random_crossings, 5);
random_p1 = prctile(random_crossings, 1);

num_as_small = sum(random_crossings <= C_sign);
share_as_small = num_as_small / num_random;

fprintf('Number of random partitions: %d\n', num_random);
fprintf('Random partitions use |L| = %d and |R| = %d.\n\n', size_L, size_R);

fprintf('Sign-split crossing weight: %.12f\n\n', C_sign);

fprintf('Random crossing weight summary:\n');
fprintf('  Mean:               %.12f\n', random_mean);
fprintf('  Standard deviation: %.12f\n', random_sd);
fprintf('  Minimum:            %.12f\n', random_min);
fprintf('  Median:             %.12f\n', random_median);
fprintf('  10th percentile:    %.12f\n', random_p10);
fprintf('  5th percentile:     %.12f\n', random_p5);
fprintf('  1st percentile:     %.12f\n', random_p1);

fprintf('\nNumber of random partitions with crossing weight <= sign split: %d\n', num_as_small);
fprintf('Share of random partitions with crossing weight <= sign split: %.6f\n\n', share_as_small);

fid = fopen(report_filename, 'a');
fprintf(fid, 'Number of random partitions: %d\n', num_random);
fprintf(fid, 'Random partitions use |L| = %d and |R| = %d.\n\n', size_L, size_R);

fprintf(fid, 'Sign-split crossing weight: %.12f\n\n', C_sign);

fprintf(fid, 'Random crossing weight summary:\n');
fprintf(fid, '  Mean:               %.12f\n', random_mean);
fprintf(fid, '  Standard deviation: %.12f\n', random_sd);
fprintf(fid, '  Minimum:            %.12f\n', random_min);
fprintf(fid, '  Median:             %.12f\n', random_median);
fprintf(fid, '  10th percentile:    %.12f\n', random_p10);
fprintf(fid, '  5th percentile:     %.12f\n', random_p5);
fprintf(fid, '  1st percentile:     %.12f\n', random_p1);

fprintf(fid, '\nNumber of random partitions with crossing weight <= sign split: %d\n', num_as_small);
fprintf(fid, 'Share of random partitions with crossing weight <= sign split: %.6f\n\n', share_as_small);
fclose(fid);


%% Final summary

fprintf('============================================================\n');
fprintf('Final Summary\n');
fprintf('============================================================\n\n');

fprintf('Part (a):\n');
fprintf('  Largest SCC size: %d out of %d original nodes.\n', largest_size, n_original);
fprintf('  Indices not in largest SCC: ');
fprintf('%d ', outside_indices);
fprintf('\n\n');

fprintf('Part (b):\n');
fprintf('  lambda_2: %.12f + %.12fi\n', real(lambda2), imag(lambda2));
fprintf('  sigma is normalized so that norm(sigma) = %.12f.\n\n', norm(sigma));

fprintf('Part (c):\n');
fprintf('  |L| = %d, |R| = %d\n', size_L, size_R);
fprintf('  L: ');
fprintf('%d ', L);
fprintf('\n');
fprintf('  R: ');
fprintf('%d ', R);
fprintf('\n\n');

fprintf('Part (d):\n');
fprintf('  C(L,R) for sign split: %.12f\n', C_sign);
fprintf('  Random mean crossing weight: %.12f\n', random_mean);
fprintf('  Random standard deviation: %.12f\n', random_sd);
fprintf('  Random minimum crossing weight: %.12f\n', random_min);
fprintf('  Random median crossing weight: %.12f\n', random_median);
fprintf('  Random 5th percentile crossing weight: %.12f\n', random_p5);
fprintf('  Number random <= sign split: %d out of %d\n', num_as_small, num_random);
fprintf('  Share random <= sign split: %.6f\n\n', share_as_small);

fprintf('Report saved to: %s\n\n', report_filename);

fid = fopen(report_filename, 'a');
fprintf(fid, '============================================================\n');
fprintf(fid, 'Final Summary\n');
fprintf(fid, '============================================================\n\n');

fprintf(fid, 'Part (a):\n');
fprintf(fid, '  Largest SCC size: %d out of %d original nodes.\n', largest_size, n_original);
fprintf(fid, '  Indices not in largest SCC: ');
fprintf(fid, '%d ', outside_indices);
fprintf(fid, '\n\n');

fprintf(fid, 'Part (b):\n');
fprintf(fid, '  lambda_2: %.12f + %.12fi\n', real(lambda2), imag(lambda2));
fprintf(fid, '  sigma is normalized so that norm(sigma) = %.12f.\n\n', norm(sigma));

fprintf(fid, 'Part (c):\n');
fprintf(fid, '  |L| = %d, |R| = %d\n', size_L, size_R);
fprintf(fid, '  L: ');
fprintf(fid, '%d ', L);
fprintf(fid, '\n');
fprintf(fid, '  R: ');
fprintf(fid, '%d ', R);
fprintf(fid, '\n\n');

fprintf(fid, 'Part (d):\n');
fprintf(fid, '  C(L,R) for sign split: %.12f\n', C_sign);
fprintf(fid, '  Random mean crossing weight: %.12f\n', random_mean);
fprintf(fid, '  Random standard deviation: %.12f\n', random_sd);
fprintf(fid, '  Random minimum crossing weight: %.12f\n', random_min);
fprintf(fid, '  Random median crossing weight: %.12f\n', random_median);
fprintf(fid, '  Random 5th percentile crossing weight: %.12f\n', random_p5);
fprintf(fid, '  Number random <= sign split: %d out of %d\n', num_as_small, num_random);
fprintf(fid, '  Share random <= sign split: %.6f\n\n', share_as_small);

fprintf(fid, 'Report saved to: %s\n\n', report_filename);
fclose(fid);


%% Optional: Plot histogram of random crossing weights

figure;
histogram(random_crossings, 40);
hold on;
xline(C_sign, 'LineWidth', 2);
hold off;

xlabel('Crossing weight');
ylabel('Frequency');
title('Random crossing weights compared with sign split');

fprintf('Histogram plotted. The vertical line is the sign-split crossing weight.\n');