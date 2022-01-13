% CCH_CONV                  predictor and p-values for CCH using convolution
%
% CALL                      [ PVALS, PRED, QVALS ] = CCH_CONV( CCH, W )
%
% GETS                      CCH         vector (a single CCH) or matrix (CCHs in columns)
%                                           has to be non-negative integers (counts)
%                           W           window width [samples] {5}
%                                           has to be non-negative integer no larger than the CCH length
% RETURNS                   PVALS       p-values (bin-wise) for exceeding chance
%                           PRED        predictor(expected values)
%                           QVALS       p-values (bin-wise) for below chance
%
% ADVICE                    for minimal run-time, collect multiple CCHs in
%                               the columns of CCH and call this routine once
%
% ADVANCED CALL             [ PVALS, PRED ] = CCH_CONV( CCH, W, WINTYPE, HF )
%
% ADVANCED ARGUMENTS        WINTYPE     window type;
%                                           {'gaussian'}    with SD of W/2; has optimal statistical properties
%                                           'rect'          of W samples; equivalent to jittering one spike train by a rectangular window of width W
%                                           'triang'        of 2W samples; equivalent to jittering both trains by a rectangular window of width W
%                                           'median'        of W samples
%                           HF          hollowed fraction; default value depends on window type;
%                                           gaussian        0.6
%                                           rectangular     0.42
%                                           triangular      0.63
%                                           median          0
%
% REFERENCE                 (1) Stark and Abeles JNM 2009
%                           (2) Spivak, Levi, Sloin, Someck, and Stark 2022

% 12-aug-09 ES

% revisions
% 11-jan-12 added qvals for deficient counts. to get global significance
%           (including correction for multiple comparisons), check crossing of alpha
%           divided by the number of bins tested
% 11-jan-22 implemented (possibly hollowed) median filter locally

function [ pvals, pred, qvals ] = cch_conv( CCH, W, WINTYPE, HF, CALCP )

% 1. CHECK ARGUMENTS
nargs                           = nargin;
if nargs < 1
    error( 'missing argument CCH' )
end
[ m, n ]                        = size( CCH );
if m * n <= 1
    error( 'improper argument CCH' )
end
if m == 1
    CCH                         = CCH'; 
    nsamps                      = n;
    ncchs                       = 1;
else
    nsamps                      = m;
    ncchs                       = n;
end
if ( sum( sum( CCH < 0 ) ) > 0 )
    error( 'improper argument CCH (must contain non-negative integers)' )
end
    
if nargs < 2 || isempty( W )
    W                           = 5; 
end
if W ~= round( W ) || W < 1
    error( 'W must be non-negative interger' )
end
if W > m
    error( 'W cannot be larger than number of CCH rows' )
end

if nargs < 3 || isempty( WINTYPE )
    WINTYPE                     = 'gaussian'; 
end
WINTYPE                         = lower( WINTYPE );

if nargs < 4 || isempty( HF )
    switch WINTYPE
        case { 'gauss', 'gaussian' }
            HF                  = 0.6;
        case 'rect'
            HF                  = 0.42;
        case 'triang'
            HF                  = 0.63;
        case 'median'
            HF                  = 1;
    end
else
    if HF < 0 || HF > 1
        error( 'HF not in range (0-1)' )
    end    
end

if nargs < 5 || isempty( CALCP )
    CALCP                       = 1;
end

% 2. PREPARE THE CONVOLUTION WINDOW
switch WINTYPE
    case { 'gauss', 'gaussian' }
        SDG                     = W / 2;
        if round( SDG ) == SDG % even W
            win                 = local_gausskernel( SDG, 6 * SDG + 1 );
            cidx                = SDG * 3 + 1;
        else
            win                 = local_gausskernel( SDG, 6 * SDG + 2 ); 
            cidx                = SDG * 3 + 1.5;
        end
    case { 'rect', 'median' }
        if W / 2 == floor( W / 2 ) % even
            win                 = ones( 1, W + 1 );
            cidx                = W / 2 + 1;
        else
            win                 = ones( 1, W );
            cidx                = ceil( W / 2 );
        end
    case 'triang'
        if W / 2 == floor( W / 2 ) % even
            win                 = triang( 2 * W + 1 );
            cidx                = W + 1;
        else
            win                 = triang( 2 * W - 1 );
            cidx                = W;
        end
    otherwise
        error( 'un-supported window type' )
end
win( cidx )                     = win( cidx ) * ( 1 - HF );
win                             = win / sum( win );      
if nsamps < ( sqrt( 2 ) * length( win ) )
    error( 'CCH-W mismatch (CCHs should be in columns; otherwise reduce W or elongate CCH)' )
end

% 3. COMPUTE A PREDICTOR BY CONVOLVING THE CCH WITH THE WINDOW:
if isequal( WINTYPE, 'median' )
    pred                        = local_medfilt( CCH, W, HF );              % hollowed non-linear filtering with edge mirroring
else
    pred                        = local_firfilt( CCH, win );                
end

% 4. COMPUTE P-VALUE BASED ON A POISSON DISTRIBUTION WITH A CONTINUITY CORRECTION:
if CALCP 
    CCH                         = round( CCH );
    pvals                       = 1 - poisscdf( CCH - 1, pred ) - poisspdf( CCH, pred ) .* ones( nsamps, ncchs ) * 0.5; % excess, deterministic
else
    pvals                       = NaN;
end
qvals                           = 1 - pvals; % deficient

return

%------------------------------------------------------------------------
% LOCAL FUNCTIONS
function Y = local_firfilt( x, W )                                          % zero-phase lag low-pass filtering of x's columns with the FIR W
C                               = length( W );
D                               = ceil( C / 2 ) - 1;
Y                               = filter( W, 1, [ flipud( x( 1 : C, : ) ); x; flipud( x( end - C + 1 : end, : ) ) ] );
Y                               = Y( 1 + C + D : end - C + D, : );
return

function K = local_gausskernel( sigmaX, N )                                 % 1D Gaussian kernel K with N samples and SD sigmaX
x                               = -( N - 1 ) / 2 : ( N - 1 ) / 2;
K                               = 1 / ( 2 * pi * sigmaX ) * exp( -( x.^2 / 2 / sigmaX^2 ) );
return

function y = local_medfilt( x, W, hf )                                      % median filtering with edge mirroring and optional hollowing
[ m, n ]                        = size( x );
hf                              = hf > 0.5;
W                               = floor( W / 2 ) * 2 + 1;
hw                              = floor( W / 2 );
x1                              = x( hw : -1 : 1, : );
x2                              = x( m : -1 : m - hw + 1, : );
x                               = [ x1; x; x2 ];
y                               = zeros( m, n );
new_sidx                        = zeros( W - 1, 1 );
new_sx                          = zeros( W - 1, 1 );
hsx                             = zeros( W - 1, 1 );
for c                           = 1 : n
    for r                       = 1 : m
        if r == 1
            xidx                = x( 1 : W, c );
            [ sx, sidx ]        = sort( xidx );
        else
            j                   = 0;
            for i               = 1 : W
                if sidx( i ) ~= 1
                    j           = j + 1;
                    new_sidx( j ) = sidx( i ) - 1;
                    new_sx( j ) = sx( i );
                end
            end
            nx                  = x( r + W - 1, c );
            j                   = 1;
            if nx > new_sx( 1 )
                while j < W
                    if nx <= new_sx( j )
                        break
                    else
                        j       = j + 1;
                    end
                end
            end
            k                   = 1;
            for i               = 1 : W
                if i == j
                    sx( i )     = nx;
                    sidx( i )   = W;
                else
                    sx( i )     = new_sx( k );
                    sidx( i )   = new_sidx( k );
                    k           = k + 1;
                end
            end
        end
        if hf
            j                   = 0;
            for i               = 1 : W
                if sidx( i ) ~= ( hw + 1 )
                    j           = j + 1;
                    hsx( j )    = sx( i );
                end
            end
            mx                  = ( hsx( hw ) + hsx( hw + 1 ) ) / 2;
        else
            mx                  = sx( hw + 1 );
        end
        y( r, c )               = mx;
    end
end

return

% EOF

