function assig = canonical_labels( labels, n )

% function assig = canonical_labels( labels, n )
%
% Maps a set of arbitrary cluster labels to labels in 1,2,...k
%
% labels( n ) = a vector of cluster labels with elements in a set of k
%		values
% assig( n ) = vector of cluster labels taking values in 1,2,... k

[ ss, isort ] = sort( labels );
[ dummy, invsort ] = sort( isort ); % the inverse permutation

cumsizes = [ find( diff( ss )) n ];
k = length( cumsizes );
assig = ones( size( labels ));
for ik = 1:k-1;
    assig(( cumsizes( ik )+1 ):cumsizes( ik+1 ) ) = ik+1;
end;
assig = assig( invsort );





