function [c, cluster_map, sigma, rho] = compress_sample(c, sigma, rho)
    cluster_map = uint32(unique(c));
    sigma = sigma(:, cluster_map + 1);
    rho = rho(:, cluster_map + 1);