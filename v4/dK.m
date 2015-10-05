function d = dK( perm1, perm2, n, nitems )

% function d = dK( perm1, perm2, n, nitems )
%
% Computes the Kendall distance between two (partial) rankings
% where the items are relabeled with arbitrary numbers.
%
% Does not work if the rankings don't contain the same items.
%
% n 	   = the maximum item appearing in perm1, perm2
% nitems   = length of perm1, perm2
% 
%
 global Q0   % dimension of Q0 varies! it has to be sufficiently large

%Q0 = triu( ones( n ), 1 );

[dum, isort1] = sort( perm1 );
[dum, isort2] = sort( perm2 );

d = nitems*(nitems-1)/2 - sum(sum(Q0(isort1,isort1).*Q0(isort2,isort2)));




