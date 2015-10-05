inv_pi = ...
    [   ...
        3   4   0   1   5
        3   1   0   2   7
        4   2   0   0   0
        2   5   7   0   0
    ]';
t = [   5	5	3	3   ]';
n = 8;

r_0 = 1 * ones(max(t), 1);
nu_0 = 1;
alpha_0 = 1;
initclusters = 2;
sigma_gibbs = 10;
rho_slices = 10;
start = 1;
iterations = 100;
thinning = 20;
filename = 'temp';
resumefilename = '';
temperatureSchedule = ones(iterations, 1);
model = 1;
randseed = 0;




sample = inference(inv_pi, t, n, r_0, nu_0, alpha_0, ...
                   initclusters, sigma_gibbs, rho_slices, iterations, thinning, ...
                   filename, resumefilename, temperatureSchedule, model, randseed);
