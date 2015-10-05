function [rho] = sample_rho(coeff, nminusj, nu, lastRho, burnin)

    function prob = logposterior(rho)
        if rho > 0
            prob = -rho * coeff - nu * (log1p(-exp(-nminusj * rho)) - log1p(-exp(-rho)));
        else
            prob = -Inf;
        end
    end

    function prob = neglogposterior(rho)
        prob = -logposterior(rho);
    end

% Start around mode
initRho = fminbnd(@neglogposterior, 0, 10, optimset('TolX', 0.1, 'Display', 'off'));
rho = slicesample(initRho, 1, 'logpdf', @logposterior, 'width', 5, 'burnin', burnin);

end