% cch_stg                   estimate STG from CCH and ACHs
% 
% call                      [ eSTG1, eSTG2, act, sil, dcCCH, crCCH ] = cch_stg( cch, ach1, nspks1, ach2, nspks2, cchbins, W )
%
% gets                      cch             [counts]    count CCHs, (nbins x npairs)
%                           ach1            [counts]    count ACHs, for the trigger (reference) unit in the CCHs
%                                                       zero-lag bin should be set to zero, (nbins x npairs)
%                           nspks1                      total number of spikes in the raw spike train from which ach1 was computed
%                           ach2            [counts]    same as ach1, for the referrred unit
%                           nspks2                      same as nspks1, for the referrred unit
%                           cchbins         [s]         nbins x 1, vector of time lags
%                           W               [samples]   span of median filter (2*delta + 1)
%
% returns                   eSTG1           STG of positive extremum in the ROI
%                           eSTG2           STG of negative extremum in the ROI
%                           act             excitation detected
%                           sil             inhibition detected
%                           dcCCH           [counts]    deconvolved CCHs
%                           crCCH           [spk/s]     conditional rate CCHs
%                           
% does 
% (1) CCH -> dcCCH:         deconvolve ACH1 and ACH2 from the CCH (if an ACH is all zeros, nothing will be deconvolved)
% (2) dcCCH -> pred:        compute the hollowed median filter predictor
% (3) dcCCH -> crCCH:       compute the conditional rate CCH
% (4) crCCH -> eSTGs:       compute the estimated STGs from the crCCH
% (5) dcCCH -> act/sil:     determine if any bin in the ROI in the dcCCH is sign. high/low
%
% note                      if CCH, ACH1, and ACH2 are matrices, 
%                               eSTG1/eSTG2/act/sil will be vectors
%                               
% calls                     cchdeconv, cch_conv, calc_stg
%
% call example              [ g1, g2, act, sil, dcCCH, crCCH, cchbins ] = cch_stg( cch, ach1, nspks1, ach2, nspks2, t * BinSizeMS / 1000, 11 );
%
% reference                 Spivak, Levi, Sloin, Someck, and Stark 2022

% 11-jan-22 ES

function [ g1, g2, act, sil, dcCCH, crCCH ] = cch_stg( cch, ach1, nspks1, ach2, nspks2, cchbins, W )

%-------------------------------------------------------------------------
% constants
% computation of the predictor (step 2)
convType                        = 'median';
hollowF                         = 1;                                        % hollowed

% detection and quantification (steps 4 and 5)
roiMS                           = [ NaN 5 ];                                % causal range until 5 ms

% detection (step 5)
alfa                            = 0.001;

%-------------------------------------------------------------------------
% argument handling
nargs                           = nargin;
if nargs < 6 || isempty( cch ) || isempty( ach1 ) || isempty( nspks1 ) || isempty( ach2 ) || isempty( nspks2 ) || isempty( cchbins )
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
if isvector( cchbins )
    cchbins                     = cchbins( : );
end
[ m,  n0 ]                      = size( cch );
[ ma1, na1 ]                 	= size( ach1 );
[ ma2, na2 ]                 	= size( ach2 );
[ mb, nb ]                      = size( cchbins );
if ( m - 1 ) / 2 ~= round( ( m - 1 ) / 2 )
    error( 'cch must have an odd number of rows (time samples)' )
end
if ~isequal( m, ma1, ma2, mb )
    error( 'input size mismatch: cch, ach1, ach2, and cchbins must have same number of rows (time samples)' )
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
if nb ~= 1
    error( 'input size mismatch: cchbins must be a vector' )
end
if nargs < 7 || isempty( W )
    W                           = 11;                                       % span of median filter; order = W - hollowF
end
W                               = W( 1 );
if ( W - 1 ) / 2 ~= round( ( W - 1 ) / 2 )
    error( 'W must be odd' )
end

%-------------------------------------------------------------------------
% preparations
dt                              = diff( cchbins( 1 : 2 ) );                 % [s]
if isnan( roiMS( 1 ) )
    roiMS( 1 )                  = ( dt * 1000 ) / 2;                        % causality imposed
end

%-------------------------------------------------------------------------
% compute

% (1) deconvolve ACHs from CCH (Eq. 21.2)
dcCCH                           = cchdeconv( cch, ach1, nspks1, ach2, nspks2 );

% (2) compute the hollowed median filter predictor (Eq. 3.3)
[ ~, pred ]                     = cch_conv( dcCCH, W, convType, hollowF, 0 );

% (3) compute conditional rate CCH (Eq. 4)
den                             = ( ones( m, 1 ) * nspks1 ) * dt; 
crCCH                           = ( dcCCH - pred ) ./ den;                  % [spk/s]

% (4) compute the STGs (Eq. 5)
t                               = cchbins * 1000;                           % [ms]
t_ROI                           = t >= roiMS( 1 ) & t <= roiMS( 2 );
[ g1, g2 ]                      = calc_stg( crCCH, t_ROI, dt );

% (5) determine if any bin in the ROI is significant
nBonf                           = sum( t_ROI );
gbUpper                         = poissinv( 1 - alfa / nBonf, max( pred( t_ROI, : ), [], 1 ) );
gbLower                         = poissinv( alfa / nBonf, min( pred( t_ROI, : ), [], 1 ) );
act                             = any( dcCCH( t_ROI, : ) >= ones( sum( t_ROI ), 1 ) * gbUpper, 1 );
sil                             = any( dcCCH( t_ROI, : ) <= ones( sum( t_ROI ), 1 ) * gbLower, 1 );

return

% EOF

