function [ pp v cdftable ] = sample_from_theta( theta, n, nsamples, cdftable, isV, t )

% function [ pp v cdftable ] = sample_from_theta( theta, n, nsamples, cdftable, isV, t )
%
% SAMPLE_FROM_THETA Takes nsamples from the distribution dictated by
% the given theta
% pp( nsamples, n )  = sample permutations, one per row
% v( n-1, nsamples ) = same permutations in the V or S representation
%		       v( j, i ) = V_j+1 in permutation i
% cdftable( n, n )   = the distributions of the V's
%		       reconstructed when empty
%
% theta( n-1 ) or (t)= parameter vector
% nsamples           = number samples required
% isV                = V or S representation
% For now, no truncation!!
%       TO DO: different sampling for large n-j
% t		     t < n truncate permutations at t (can't be missing)
%		     = n full permutations


% Create cdftable

if isempty( cdftable )  % construct probability table
   hh = hankel(ones(n,1));
   table = exp(-[theta 1]'*(0:(n-1)));  
   if ~isV
      if t < n
          hh = hh( 1:t, : );
          table = table( 1:t, : );
      end;
   end;
   table = table.*hh;

   table = table./repmat(sum(table, 2), [ 1, n ]);
   cdftable = cumsum( table, 2 ).*hh;
end;  

if isV 
   pp = zeros( nsamples, n );    
   v = zeros(nsamples, n);
else
   pp = zeros( nsamples, t );
   v = zeros(nsamples, t);
end;

%  Sample v 
  
if isV   % V representation
    for j=1:n-1 
        r = rand( nsamples, 1 );
        v( :, j ) = sum( (repmat( r, [1, n-j] ) - repmat( cdftable( j, 1:n-j), [nsamples, 1]) > 0 ), 2 );
    end;
else    % S representation
    for j=1:min( t, n-1 ); 
        r = rand( nsamples, 1 );
        v( :, j ) = sum( (repmat( r, [1, n-j] ) - repmat( cdftable( j, 1:n-j), [nsamples, 1]) > 0 ), 2 );
    end;
end;


% Construct the actual ranking in the standard vector representation

v = v+1;
if isV   % v is the V representation
   iii = repmat( 1:n, [ nsamples, 1 ]);
   for j = 1:n;    
       vj = v( :, j );
       for isam = 1:nsamples;
           pp( isam, iii( isam, vj( isam ))) = j;
           iii( isam, vj( isam ):end-1) = iii( isam, vj( isam )+1:end);
       end;
       iii = iii( :, 1:end-1 );
   end;
   if  t< n  % Truncate to t
       pp = pp( :, 1:t );
   end;
   
else    % v is the S representation
    if t< n  % Truncate to t
        iii = repmat( (1:n)', [ 1, nsamples ]);  % TAKES A LOT OF MEMORY
        for j = 1:t;
            sj = v( :, j );
            sj = sj + (0:nsamples-1)'*(n-j+1); % turn them into indices into iii
            pp( :, j ) = iii( sj );
            iii = iii( setdiff( 1:nsamples*(n-j+1), sj )); % delete the occupied
                                              % position from iii
        end;
    else
        % full permutations
        iii = repmat( (1:n)', [ 1, nsamples ]);  
        for j = 1:n-1;
            sj = v( :, j );
            sj = sj + (0:nsamples-1)'*(n-j+1); % turn them into indices into iii
            pp( :, j ) = iii( sj );
            iii = iii( setdiff( 1:nsamples*(n-j+1), sj )); % delete the occupied
                                              % position from iii
        end;
        pp( :, n ) = iii;
    end;   
end;



    


   
