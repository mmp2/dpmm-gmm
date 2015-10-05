function run_cao(filename, alpha_0, initclusters)


% Read in course data, run sampler

load('../../data/course-preferences-brendan-murphy/caointegers.mat');

% Transform output to format for Harr's sampler
inv_pi = cao' - 1;
inv_pi(inv_pi < 0) = 0;
t = tt;
r_0 = 1 * ones(max(t), 1);
nu_0 = 1;

sigma_gibbs = 10;
rho_slices = 10;
iterations = 10000;
thinning = 1;
filename = strcat(filename, '_alpha', int2str(alpha_0), '_initclusters', int2str(initclusters));
resumefilename = '';
temperatureSchedule = ones(iterations, 1);
model = 0;
randseed = 0;

sample = inference(inv_pi, t, n, r_0, nu_0, alpha_0, ...
                   initclusters, sigma_gibbs, rho_slices, iterations, thinning, ...
                   filename, resumefilename, temperatureSchedule, model, randseed);

% Plot results
ctrace = zeros(N, iterations + 1);
for i = 1:(iterations + 1)
    ctrace(:, i) = canonical_labels(double(sample(i).c') + 1, N);
end

[ cluss, isort ]=sort( cluststar );

fontsize = 16;
imagesc( [ cluss; ctrace( isort, : )'; cluss ] );
imagesc(  mod([ cluss; ctrace( isort, : )'; cluss ], 8) );
set(gca, 'FontSize', fontsize );
ht = title( 'cluster labels' );
set(ht, 'FontSize', fontsize );
ht = ylabel( 'iterations' );
set(ht, 'FontSize', fontsize );
