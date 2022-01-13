% calc_stg          spike transmission gain from conditional rate CCH
% 
% call              [ g1, g2 ] = calc_stg( crcch, t_ROI, dt, stgMode, roiHardness )
%
% gets              crcch           [spks/s]; may be a matrix (nbins x npairs)
%                   t_ROI           logical vector (nbins x 1)
%                   dt              [s], CCH bin size
%                   stgMode         {1}         0   detects global maximum in ROI (always exists, may be non-causal)
%                                               1   detects local maximum in ROI (may not exist, but preserves causality)
%                   roiHardness     {[1 0]}     hard/soft edges (limits/does not limit STG base by ROI)
% 
% returns           g1              STG of positive extremum in the ROI
%                   g2              STG of negative extremum in the ROI
%
% note              if crcch is a matrix, g1 and g2 will be row vectors
% 
% calls             nothing

% 07-mar-21 ES

% revisions
% 11-jan-22 support matrix input

function [ g1, g2 ] = calc_stg( crcch, t_ROI, dt, stgMode, roiHardness )

% argument handling
nargs                           = nargin;
if nargs < 4 || isempty( stgMode )
    stgMode                     = 1;
end
if nargs < 5 || isempty( roiHardness )
    roiHardness                 = [ 1 0 ];
end

% prepare
ft_ROI                          = find( t_ROI );
n                               = size( crcch, 1 );
crcch( isnan( crcch ) )           = 0;

% check input size
[ m1, n1 ]                      = size( crcch );
[ m2, n2 ]                      = size( t_ROI );
if m1 ~= m2
    error( 'input size mismatch: crcch and t_ROI must have same number of rows' )
end
if n2 ~= 1
    error( 'input size error: t_ROI must have a single column' )
end
t_ROI                           = t_ROI * ones( 1, n1 );

% run recursively if matrices
if n1 > 1
    g1                          = NaN( 1, n1 );
    g2                          = g1;
    for c = 1 : n1
        [ g1( c ), g2( c ) ]     = calc_stg( crcch( :, c ), t_ROI( :, c ) ...
            , dt, stgMode, roiHardness );
    end
    return
end
    
% support of positive extremum
switch stgMode
    case 0
        x                       = crcch( t_ROI );
        [ ~, maxidx ]           = max( x );
    case 1
        pidx                    = ft_ROI( [ 1 end ] ) + [ -1 1 ]';
        xidx                    = [ pidx( 1 ); ft_ROI; pidx( 2 ) ];
        nROI                    = length( ft_ROI ) + 1;
        x                       = crcch( xidx );
        sx                      = find( diff( x ) == 0 );                   % same values of x
        ux                      = x;
        ix                      = ( 1 : length( x ) )';
        ix( sx )                = [];
        ux( sx )                = [];
        tmpidx                  = find( diff( sign( diff( ux ) ) ) < -1 ) + 1;
        maxidx                  = ix( tmpidx ) - 1;
        if ~isempty( sx ) && sx( end ) == nROI && length( ux ) > 1 && ux( end ) > ux( end - 1 )
            maxidx              = sx( end ) - 1;
        end
end
if length( maxidx ) > 1
    [ ~, subidx ]               = max( x( maxidx ) );
    maxidx                      = maxidx( subidx );
end
if isempty( maxidx )
    sidx                        = [];
else
    pidx                    	= ft_ROI( maxidx );
    si                        	= 1;                                    	% first preceding negative value
    for i                     	= pidx : -1 : 1
        if crcch( i ) < 0
            si                 	= i + 1;
            break
        end
    end
    ei                         	= n;                                    	% first proceeding negative value
    for i                      	= pidx : n
        if crcch( i ) < 0
            ei                 	= i - 1;
            break
        end
    end
    sidx                       	= si : ei;
end
if roiHardness( 1 )
    rmv                         = sidx < ft_ROI( 1 );
    sidx( rmv )                 = [];
end
if roiHardness( 2 )
    rmv                         = sidx > ft_ROI( end );
    sidx( rmv )                 = [];
end
% compute the integral
a1                              = nanmean( crcch( sidx ) );
b1                              = sum( ~isnan( crcch( sidx ) ) ) * dt;
g1                              = a1 * b1;

% support of negative extremum
switch stgMode
    case 0
        [ ~, minidx ]           = min( x );
    case 1
        tmpidx                  = find( diff( sign( diff( ux ) ) ) > 1 ) + 1;
        minidx                  = ix( tmpidx ) - 1;
        if ~isempty( sx ) && sx( end ) == nROI && length( ux ) > 1 && ux( end ) < ux( end - 1 )
            minidx              = sx( end ) - 1;
        end
end
if length( minidx ) > 1
    [ ~, subidx ]               = min( x( minidx ) );
    minidx                      = minidx( subidx );
end
if isempty( minidx )
    sidx                        = [];
else
    pidx                     	= ft_ROI( minidx );
    si                         	= 1;                                       	% first preceding positive value
    for i                      	= pidx : -1 : 1
        if crcch( i ) > 0
            si                	= i + 1;
            break
        end
    end
    ei                       	= n;                                    	% first proceeding positive value
    for i                      	= pidx : n
        if crcch( i ) > 0
            ei               	= i - 1;
            break
        end
    end
    sidx                     	= si : ei;
end
if roiHardness( 1 )
    rmv                         = sidx < ft_ROI( 1 );
    sidx( rmv )                 = [];
end
if roiHardness( 2 )
    rmv                         = sidx > ft_ROI( end );
    sidx( rmv )                 = [];
end
% compute the integral
a2                              = nanmean( crcch( sidx ) );
b2                              = sum( ~isnan( crcch( sidx ) ) ) * dt;
g2                              = a2 * b2;

return

% EOF