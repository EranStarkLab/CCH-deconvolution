% cchdeconv         deconvolve ACHs from CCHs, yielding dcCCHs
%
% call              dccch = cchdeconv( cch, ach1, nspks1, ach2, nspks2 )
%
% gets              cch                 CCH (vector or matrix), [counts]
%                   ach1                ACH of the trigger spike train, [counts]
%                   nspks1              number of spikes in the trigger train
%                   ach2                ACH of the referred train, [counts]
%                   nspks2              number of spikes in the referred train
%
% returns           dccch              	the dcCCH, [counts] (same size as CCH)
%
% does              (1) clip and scale ACHs
%                   (2) deconvolve ACHs from CCH
%
% notes
% (1) cch, ach1, and ach2 must be same-sized matrices (nbins x npairs)
% (2) nspks1 and nspks2 must have either a size of 1 of npairs
% (3) the number of time samples in the CCH, nbins, must be odd
%
% calls             nothing
%
% reference         Spivak, Levi, Sloin, Someck, and Stark 2022

% 14-oct-20 ES

% revisions
% 11-jan-22 simplified

function dccch                  = cchdeconv( cch, ach1, nspks1, ach2, nspks2 )

% argument handling
nargs                           = nargin;
if nargs < 5 || isempty( cch ) || isempty( ach1 ) || isempty( nspks1 ) || isempty( ach2 ) || isempty( nspks2 )
    error( 'missing arguments' )
end
if isvector( cch )
    cch                         = cch( : );
end
if isvector( ach1 )
    ach1                        = ach1( : );
end
if isvector( nspks1 )
    nspks1                      = nspks1( : )';
end
if isvector( ach2 )
    ach2                        = ach2( : );
end
if isvector( nspks2 )
    nspks2                      = nspks2( : )';
end
[ m,  n0 ]                      = size( cch );
[ ma1, na1 ]                 	= size( ach1 );
[ ma2, na2 ]                 	= size( ach2 );
if ( m - 1 ) / 2 ~= round( ( m - 1 ) / 2 )
    error( 'cch must have an odd number of rows (time samples)' )
end
if ~isequal( m, ma1, ma2 )
    error( 'input size mismatch: cch, ach1, and ach2 must have same number of rows (time samples)' )
end
if ~isequal( n0, na1, na2 )
    error( 'input size mismatch: cch, ach1, and ach2 must have same number of columns' )
end
if numel( nspks1 ) == 1 && n0 > 1
    nspks1                      = nspks1 * ones( 1, n0 );
end
if numel( nspks2 ) == 1 && n0 > 1
    nspks2                      = nspks2 * ones( 1, n0 );
end
if size( nspks1, 1 ) ~= 1 || size( nspks1, 2 ) ~= n0
    error( 'input size mismatch: nspks1 must have the same number of elements as the number of columns in cch' )
end
if size( nspks2, 1 ) ~= 1 || size( nspks2, 2 ) ~= n0
    error( 'input size mismatch: nspks2 must have the same number of elements as the number of columns in cch' )
end

% preparations
hw                              = ( m - 1 ) / 2;

% scale ACHs 
ach1k                           = ach1;
ach1k                           = ach1k - sum( ach1k ) / m;                 % remove mean of clipped
ach1k                           = ach1k ./ nspks1;                          % divide by nspks1
hidx                            = [ 1 : hw ( hw + 2 ) : m ];
ach1k( hw + 1, : )              = 1 - sum( ach1k( hidx, : ) );              % set zero-lag bin s.t. sum of 1

ach2k                           = ach2;
ach2k                           = ach2k - sum( ach2k ) / m;                 % remove mean of clipped
ach2k                           = ach2k ./ nspks2;                          % divide by nspks1
hidx                            = [ 1 : hw ( hw + 2 ) : m ];
ach2k( hw + 1, : )              = 1 - sum( ach2k( hidx, : ) );              % set zero-lag bin s.t. sum of 1

% deconvolve ACHs from the CCH
den                             = fft( ach1k, m ) .* fft( ach2k, m );
dccch                          	= real( ifft( fft( cch, m ) ./ den, m ) );

% organize output
dccch                           = dccch( [ 2 : m 1 ], : );                  % shift DC to end
dccch( dccch < 0 )              = 0;                                        % clip negatives (may occur numerically) to zero

return

% EOF

