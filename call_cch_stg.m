% call_cch_stg          Wrapper for cch_stg, estimate STG from spike times
%  
% call                	[ eSTG1, eSTG2, act, sil, dcCCH, crCCH,cchbins ] = call_cch_stg( spkT, spkL, spkFS, binSize, halfSize, W )
% 
% gets                	spkT                        [samples]   spike times [samples] of all units (Nspikes x 1 )
%                   	spkL                                    spikes labels of all units (Nspikes x 1 ) 
%                   	spkFS                       [samples/s] recording sampling rate
% 
%                       binSize         {0.001};    [s]         bin size for creating the dcCCH
%                       halfSize        {0.05};     [s]         half size of the dcCCH,
%                     	W               {11};       [samples]   span of median filter (2*delta + 1)
% 
% returns            	eSTG1                       [spk/spk]   STG of positive extremum in the ROI
%                    	eSTG2                       [spk/spk]   STG of negative extremum in the ROI
%                   	act                                     excitation detected
%                       sil                                     inhibition detected
%                       dcCCH                       [counts]    deconvolved CCHs
%                       crCCH                       [spk/s]     conditional rate CCHs
%                   	cchbins                     [s]         nbins x 1, vector of time lags
% does 
% (1) calls CCG 
% (2) computes Nspikes
% (3) organize all variables
% (4) calls cch_stg
%                              
% calls                 (1) CCG (FMAToolbox)
%                       (2) cch_stg (Spivak et al., 2022)
%
% reference             Spivak, Levi, Sloin, Someck, and Stark 2022

% 12-jan-22 LS

function [ eSTG1, eSTG2, act, sil, dcCCH, crCCH,cchbins ] = call_cch_stg( spkT, spkL, spkFS, binSize, halfSize, W )


%-------------------------------------------------------------------------
% argument handling
nargs                           = nargin;
if nargs < 3 || isempty( spkT ) || isempty( spkL ) || isempty( spkFS ) 
    error( 'missing arguments' )
end

if isvector( spkT )
    spkT                    	= spkT( : );
else 
    error( 'Spike times must be a vector' )
end

if isvector( spkL )
    spkL                    	= spkL( : );
else 
    error( 'Spike labels must be a vector' )
end

if length( spkL ) ~= length( spkT ) 
    error( 'Spike labels and Spike times must have the same size' )
end

if nargs < 4 || isempty( binSize )
    binSize                     = 0.001;                                    
end

if nargs < 5 || isempty( halfSize )
    halfSize                    = 0.05;                      
end

if nargs < 6 || isempty( W )
    W                           = 11;                         
end
W                               = W( 1 );
if ( W - 1 ) / 2 ~= round( ( W - 1 ) / 2 )
    error( 'W must be odd' )
end
%-------------------------------------------------------------------------------

% compute all CCGs
[ cc, cchbins ]                	= CCG( spkT./spkFS, spkL, 'binSize', binSize, 'duration', halfSize, 'norm', 'counts');

% compute Nspikes for each unit
Gsub                            = unique( spkL );
Nspk                            = zeros( size( Gsub ) );
nGroups                         = length( Gsub );
for g                           = 1 : nGroups
    Nspk( g )                   = sum( spkL == Gsub( g ) );
end

% reshape CCHs ACHs and Nspikes to be 2D matrices
cc                              = cc( :, Gsub, Gsub );
[ cch, ach1, nspks1 ...
    , ach2, nspks2 ]            = CCH3D_reshape( cc, Nspk );

% call cch_stg
[ eSTG1, eSTG2, act, sil ...
    , dcCCH, crCCH ]            = cch_stg( cch, ach1, nspks1, ach2, nspks2, cchbins, W );

return

%------------------------------------------------------------------------
% LOCAL FUNCTIONS
function [ cch, ach1, nspks1, ach2, nspks2 ]  = CCH3D_reshape( ccg, nspks0 )            % convert 3D matrix to three 2D matrices 

[ m, n, n2 ]                    = size( ccg );
if n ~= n2 
    error( 'must be a 3D matrix with equal number of colums and sheets' )
end

% reshape cch and extract ach to prepare for deconvolution
cch                             = reshape( ccg, [ m n * n ] );            % columns are [ 11 12 ... 1n 21 22 ... 2n ... ]
    
% source column indices for ACH (n x n matrix organized as an 1 x n^2 vector):
aidx                            = ( 0 : ( n - 1 ) ) * n + ( 1 : n );
aidx                            = ones( n, 1 ) * aidx;
aidx                            = aidx( : );

% target column indices for ACH1:
aidx1                           = reshape( ( 1 : n^2 ), [ n n ] )';
aidx1                           = aidx1( : );
ach1                            = zeros( m, n * n );
ach1( :, aidx1 )                = cch( :, aidx );

% target column indices for ACH2:
aidx2                           = ( 1 : n^2 )';
ach2                            = zeros( m, n * n );
ach2( :, aidx2 )                = cch( :, aidx );

% spike count indices
nidx                            = ( 1 : n )' * ones( 1, n );
nidx                            = nidx( : );
nspks1                          = zeros( 1, n * n );
nspks2                          = zeros( 1, n * n );
nspks1( aidx2 )                 = nspks0( nidx );
nspks2( aidx1 )                 = nspks0( nidx );

% prune auto-CCHs
ridx                            = 1 : ( n + 1 ) : n^2;
cch( :, ridx )                	= [];
ach1( :, ridx )                 = [];
ach2( :, ridx )                 = [];
nspks1( :, ridx )               = [];
nspks2( :, ridx )               = [];

return

% EOF
