function [ sigma, tsig ] = sample_from_postG( t, n, Rlik, trunc )

%function [ sigma, tsig ] = sample_from_postG( t, n, Rlik, trunc )
%
% Samples from the [conjugate] posterior  of sigma given theta
%
% t = unused
% theta(t)   = parameter theta. We sample from P( sigma | theta, ...)
%	       incorporated in Rlik
% n = total number items
% Rlik(n,n) = sufficient statistics  Rlik = sum_j theta_j*Rlik( :,:,j)
% trunc	    = 1 truncate when the remainder of Rlik becomes 0
%		equivalent to returning a permutation restricted to the 
%		observed items	    	      	      
%	    = 0 returns a whole permutation
% OUTPUT sigma = sampled permutation
%	 tsig  = length of significant sigma
% For simplicity, and compatibility with Lsigma, sigma is always filled at random.
% tsig is the rank of the last "observed" item in Rlik.
% for unobserved items Rlik has: rows = 0, columns = q ==> uniform distribution
% when only unobserved items are left, sample_from_postG samples a random permutation


sigrange = 1:n;
sigma = zeros( 1, n );
u = rand( 1, n-1 );  % the random numbers

for i=1:n-1;
    probsig = exp( -sum( Rlik( sigrange, sigrange), 1 ));
    Zsig = sum( probsig );
    if abs(Zsig - n+i-1) < eps;  % remainder is uniform distribution
       i=i-1;
       break;
    end;
    cumsig = cumsum( probsig );
    ii = sum(( u( i )*Zsig > cumsig ))+1;
    sigma( i ) = sigrange( ii );
    sigrange = sigrange( [ 1:ii-1 ii+1:end] );
end;

if trunc
   tsig = i;
else
   tsig = n;
end;
sigma( i+1:n ) = sigrange( randperm( n-i ));





