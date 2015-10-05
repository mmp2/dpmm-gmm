function [sample] = inference(inv_pi, ...
                              t, ...
                              n, ...
                              r_0, ...
                              nu_0, ...
                              alpha_0, ...
                              initclusters, ...
                              sigma_gibbs, ...
                              rho_slices, ...
                              iterations, ...
                              thinning, ...
                              filename, ...
                              resumefilename, ...
                              temperatureSchedule, ...
                              model, ...
                              randseed)
% Runs the DPMM for GMM sampler
%
% Inputs:
% inv_pi:
%   Top-t orderings, (max_t, N) matrix; all values are zero-based
% t:
%   Length of each permutation; (N) vector
% r_0:
%   Prior values for r; (max_t) vector
% nu_0:
%   Prior strength; scalar
% alpha_0:
%   Prior concentration parameter; scalar
% initclusters: 
%   Guaranteed exact number of clusters to initialize sampler at; set to
%   one for everything in one cluster, set to N for everything in its own
%   cluster
% sigma_gibbs:
%   Number of iterations to perform for sampling a new sigma for clusters
%   larger than one
% rho_slices:
%   When model type is 1 (slice-gibbs), number of slice sampling iterations
%   to run for sampling each rho
% iterations:
%   Number of iterations to perform
% thinning:
%   If filename is specified, intervals at which to dump file
% filename:
%   If not an empty string, saves intermediate progress to a file starting
%   with this name
% resumefilename:
%   NOT IMPLEMENTED
% temperatureSchedule:
%   NOT IMPLEMENTED - SET TO ones(iterations, 1)
% model:
%   0 for BETA-GIBBS, 1 for SLICE-GIBBS
% randseed:
%   0 for new seed, or specify a value to reproduce random values
%
% Returns:
% sample:
%   Struct array of samples at each iteration; first sample is starting
%   state
% sample(i).c:
%   Cluster assignments at iteration i - 1; uint32(N) vector
% sample(i).cluster_map:
%   cluster_map(j) = k means that column j in sigma and rho correspond to
%   cluster k; uint32(x) vector, where x is the number of valid clusters
% sample(i).sigma:
%   Centroid permutations, in RANK format (not ordering format); column k
%   corresponds to cluster_map(k); uint32(n, x) matrix
% sample(i).rho:
%   Dispersion parameters; column k corresponds to cluster_map(k);
%   double(max_t, x) matrix

% Track start time ========================================================
starttime = datestr(now); %#ok<NASGU>

% Seed random =============================================================
if randseed == 0
    randseed = sum(100 * clock);
end
rand('state', randseed); %#ok<RAND>

% Convert things to uint32 and column vectors =============================
inv_pi = uint32(inv_pi);
t = reshape(uint32(t), length(t), 1);
n = uint32(n);
sigma_gibbs = uint32(sigma_gibbs);
rho_slices = uint32(rho_slices);
r_0 = reshape(r_0, length(r_0), 1);
model = uint32(model);

% Basic corpus stats ======================================================
N = size(t, 1);
max_t = max(t);

if strcmp(resumefilename, '')
    % Initialize structure and estimates ==================================
    sample = struct('c',           cell(iterations + 1, 1), ...
                    'cluster_map', cell(iterations + 1, 1), ...
                    'sigma',       cell(iterations + 1, 1), ...
                    'rho',         cell(iterations + 1, 1));
    
    % Init c
    % Initialization strategy determines how to do initial clustering
    c = zeros(N, 1, 'uint32');
    
    % Take top clusters from a random permutation to make sure each
    % required cluster has at least one entry
    dataperm = randperm(N);
    c(dataperm(1:initclusters)) = uint32(1:initclusters) - 1;
    % Allocate rest randomly
    c(dataperm(initclusters + 1:N)) = uint32(floor(unifrnd(0, initclusters, N - initclusters, 1)));
    
    % Init sigma to random permutations
    sigma = zeros(n, N, 'uint32');
    for i = 1:N
        sigma(:, i) = uint32(randperm(n) - 1);
    end

    % Init rho to zero (its init value is never used) for BETA_GIBBS
    % For SLICE_GIBBS, we need to seed to a reasonable default, which we do
    % by drawing from the approximate BETA_GIBBS prior
    if model == 0
        rho = zeros(max_t, N);
    else
        rho = -log(betarnd(repmat(nu_0 .* r_0, 1, N), nu_0 + 1));
    end
    
    % Precompute pi_R and the beta_table for N = 1 sigma sampling
    pi_R = compute_pi_R(inv_pi, t, n);
    beta_table = compute_beta_table(nu_0, r_0, n, max_t);
        
    % MEX function takes care of initializing these caches
    nc = [];
    S = [];
    
    % Save initialization state
    [sample(1).c sample(1).cluster_map sample(1).sigma sample(1).rho] = ...
        compress_sample(c, sigma, rho);
else
    % Read parameters from MAT file
    % HACK HACK - this is broken right now do not use
    error('Resume is not properly supported right now');
    load(resumefilename); %#ok<UNRCH>
end

% Actual sampling procedure ===============================================
if thinning == 0
   thinning = iterations;
end
iter = 0;
likelihoods = zeros(iterations, 1);
while iter < iterations
    [c ...
     sigma ...
     rho ...
     nc ...
     S ...
     likelihood] ...
       = sample_model( ...
            inv_pi, ...
            t, ...
            n, ...
            c, ...
            sigma, ...
            rho, ...
            nc, ...
            S, ...
            r_0, ...
            nu_0, ...
            alpha_0, ...
            pi_R, ...
            beta_table, ...
            temperatureSchedule(iter + 1), ...
            sigma_gibbs, ...
            rho_slices, ...
            1, ...
            iter, ...
            model);
    iter = iter + 1;
    likelihoods(iter) = likelihood;
    
    [sample(iter + 1).c sample(iter + 1).cluster_map sample(iter + 1).sigma sample(iter + 1).rho] = ...
        compress_sample(c, sigma, rho);

    if strcmp(filename, '') == 0
        if mod(iter, thinning) == 0
            savetime = datestr(now); %#ok<NASGU>
            fprintf('Saving sample at iteration %f\n', iter);
            save(strcat(filename, '.temp.mat'), 'sample', 'nc', 'S');
            
            % If we reached here that means we probably saved successfully,
            % so overwrite old file now
            if exist(strcat(filename, '.mat'), 'file')
                delete(strcat(filename, '.mat'));
            end
            movefile(strcat(filename, '.temp.mat'), strcat(filename, '.mat'));
        end
    end
end

if strcmp(filename, '') == 0
    savetime = datestr(now); %#ok<NASGU>
    fprintf('Saving final sample\n');
    save(strcat(filename, '.mat'), 'sample', 'nc', 'S');
end