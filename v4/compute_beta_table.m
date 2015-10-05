function B = compute_beta_table(nu_0, r_0, n, t)
% Produces the probability table used for sampling sigma for N = 1.
% nu_0:     prior strength, double scalar
% r_0:      prior values, double(t) vector
% n:        number of items
% t:        maximum rank for consideration
% Output:   double(n, t) matrix matrix

    disp('Computing Beta table');
    tic;
    k = repmat((0:(double(n) - 1))', 1, t);
    r_0 = reshape(r_0, 1, t);
    B = beta(nu_0 * repmat(r_0, n, 1) + k, nu_0 + 2);
    toc;
    disp('Done: Computing Beta table');
end

