% cch_stg_demo                  load example simulated and real data and plot CCHs and dcCCHs
%
% requires:
% routines                      call_cch_stg, cch_stg, cchdeconv, cch_conv, calc_stg
% datasets                      simData.mat, CA1_Data.mat
% external sources              CCG (and dependencies; see FMAToolbox)
%
% 13-jan-22 LS

function cch_stg_demo

load( 'simData.mat', 'simData' )
load( 'CA1_Data.mat', 'CA1_Data' )

%-------------------------------------------------------------------------------------
% 0. Global parameters
%-------------------------------------------------------------------------------------
% parameters for CCH computations
binSize                         = 0.001;                                    % [s]
halfSize                        = 0.1;                                      % [s]
W                               = 11;                                       % [samples]

% constants for STG computations
roiMS                           = [ NaN 5 ]; 
convType                        = 'median';
hollowF                         = 1; 


%-------------------------------------------------------------------------------------
% 1. Simulated data
%-------------------------------------------------------------------------------------
Cbars                           = [161, 136, 127]./255;
Cdc                             = [ 0.18   0.49     0.196];
figure('Renderer', 'painters', 'Position', [10 10 900 600])


%-------------------------------------------------------------------------------------
% 1.1. focus on a simulated pair with 
% -strong mono-excitatory connection
% -high burst spiking activity

spkT                            = simData( 1 ).T;
spkL                            = simData( 1 ).label;
spkFS                           = simData( 1 ).spkFs;
rSTG                            = simData( 1 ).rSTG;

% 1.1.1. compute eSTGs from dcCCH
[ eSTG1, ~, ~, ~, dcCCH, ~, cchbins ] = call_cch_stg( spkT, spkL, spkFS, binSize, halfSize, W );

% 1.1.2. compute CCG
cc                              = CCG( spkT ./ spkFS, spkL, 'binSize', binSize, 'duration', halfSize );
Gsub                            = unique( spkL );
cc                              = cc( :, Gsub, Gsub );
bs                              = length( Gsub ) - 1;
u1                              = 1;
u2                              = 2;

% 1.1.3. compute and plot eSTGs from CCH
nspks1                          = sum( spkL == Gsub( u1 ) );
dt                              = diff( cchbins( 1 : 2 ) );                 % [s]
if isnan( roiMS( 1 ) )
    roiMS( 1 )                  = ( dt * 1000 ) / 2;                        % causality imposed
end
% (1) compute the hollowed median filter predictor (Eq. 3.3)
[ ~, pred ]                     = cch_conv( cc( :, u1, u2 ), W, convType, hollowF, 0 );

% (2) compute conditional rate CCH (Eq. 4)
den                             =  nspks1  * dt; 
crCCH                           = ( cc( :, u1, u2 ) - pred ) ./ den;      	% [spk/s]

% (3) compute the STGs (Eq. 5)
t                               = cchbins * 1000;                           % [ms]
t_ROI                           = t >= roiMS( 1 ) & t <= roiMS( 2 );
g1                              = calc_stg( crCCH, t_ROI, dt );

% (4) plot CCH and dcCCH for the selected pair
subplot( 2, 2, 1 );
plot_CCH( cchbins, cc( :, u1, u2 ), Cbars, rSTG, g1 )

subplot( 2, 2, 2 );
cidx                            = ( u2 - 1 ) * bs + u1;
plot_CCH( cchbins, dcCCH( :, cidx ), Cdc, rSTG, eSTG1( cidx ) )

%-------------------------------------------------------------------------------------
% 1.2. focus on a simulated pair with 
% -weak mono-excitatory connection
% -high burst spiking activity

spkT                            = simData( 2 ).T;
spkL                            = simData( 2 ).label;
spkFS                           = simData( 2 ).spkFs;
rSTG                            = simData( 2 ).rSTG;

% 1.2.1. compute eSTGs from dcCCH
[ eSTG1, ~, ~, ~, dcCCH, ~, cchbins ] = call_cch_stg( spkT, spkL, spkFS, binSize, halfSize, W );

% 1.2.2. compute CCG
cc                              = CCG( spkT ./ spkFS, spkL, 'binSize', binSize, 'duration', halfSize );
Gsub                            = unique( spkL );
cc                              = cc( :, Gsub, Gsub );
bs                              = length( Gsub ) - 1;
u1                              = 1;
u2                              = 2;

% 1.2.3. compute and plot eSTGs from CCH
nspks1                          = sum( spkL == Gsub( u1 ) );
dt                              = diff( cchbins( 1 : 2 ) );                 % [s]
if isnan( roiMS( 1 ) )
    roiMS( 1 )                  = ( dt * 1000 ) / 2;                        % causality imposed
end
% (1) compute the hollowed median filter predictor (Eq. 3.3)
[ ~, pred ]                     = cch_conv( cc( :, u1, u2 ), W, convType, hollowF, 0 );

% (2) compute conditional rate CCH (Eq. 4)
den                             =  nspks1  * dt; 
crCCH                           = ( cc( :, u1, u2 ) - pred ) ./ den;      	% [spk/s]

% (3) compute the STGs (Eq. 5)
t                               = cchbins * 1000;                           % [ms]
t_ROI                           = t >= roiMS( 1 ) & t <= roiMS( 2 );
g1                              = calc_stg( crCCH, t_ROI, dt );

% (4) plot CCH and dcCCH for the selected pair
subplot( 2, 2, 3 );
plot_CCH( cchbins, cc( :, u1, u2 ), Cbars, rSTG, g1 )

subplot( 2, 2, 4 );
cidx                            = ( u2 - 1 ) * bs + u1;
plot_CCH( cchbins, dcCCH( :, cidx ), Cdc, rSTG, eSTG1( cidx ) )

%-------------------------------------------------------------------------------------
% 2. Real data CA1
%-------------------------------------------------------------------------------------
Cbars                           = [161, 136, 127]./255;
Cdc                             = [ 0.18   0.49     0.196];
figure('Renderer', 'painters', 'Position', [10 10 900 600])
spkT                            = CA1_Data( 1 ).T;
spkL                            = CA1_Data( 1 ).label;
spkFS                           = CA1_Data( 1 ).spkFs;

%-------------------------------------------------------------------------------------
% 2.1. compute eSTGs from dcCCH for all pairs
[ eSTG1, ~, ~, ~, dcCCH, ~, cchbins ] = call_cch_stg( spkT, spkL, spkFS, binSize, halfSize, W );

% 2.2. compute CCGs for all pairs
cc                              = CCG( spkT ./ spkFS, spkL, 'binSize', binSize, 'duration', halfSize );
Gsub                            = unique( spkL );
cc                              = cc( :, Gsub, Gsub );
bs                              = length( Gsub ) - 1;

% 2.3. compute and plot eSTGs from CCHs for selected pairs
%-------------------------------------------------------------------------------------
% 2.3.1. compute eSTGs from CCHs for a selected pair
u1                              = 4;
u2                              = 6;

nspks1                          = sum( spkL == Gsub( u1 ) );
dt                              = diff( cchbins( 1 : 2 ) );                 % [s]
if isnan( roiMS( 1 ) )
    roiMS( 1 )                  = ( dt * 1000 ) / 2;                        % causality imposed
end
% (1) compute the hollowed median filter predictor (Eq. 3.3)
[ ~, pred ]                     = cch_conv( cc( :, u1, u2 ), W, convType, hollowF, 0 );

% (2) compute conditional rate CCH (Eq. 4)
den                             =  nspks1  * dt; 
crCCH                           = ( cc( :, u1, u2 ) - pred ) ./ den;      	% [spk/s]

% (3) compute the STGs (Eq. 5)
t                               = cchbins * 1000;                           % [ms]
t_ROI                           = t >= roiMS( 1 ) & t <= roiMS( 2 );
g1                              = calc_stg( crCCH, t_ROI, dt );

% (4) plot CCH and dcCCH for the selected pair
subplot( 2, 2, 1 );
plot_CCH( cchbins, cc( :, u1, u2 ), Cbars, NaN, g1 )

subplot( 2, 2, 2 );
cidx                            = ( u2 - 1 ) * bs + u1;
plot_CCH( cchbins, dcCCH( :, cidx ), Cdc, NaN, eSTG1( cidx ) )


%-------------------------------------------------------------------------------------
% 2.3.2. compute eSTGs from CCHs for a second pair
u1                              = 3;
u2                              = 5;

nspks1                          = sum( spkL == Gsub( u1 ) );
dt                              = diff( cchbins( 1 : 2 ) );                 % [s]
if isnan( roiMS( 1 ) )
    roiMS( 1 )                  = ( dt * 1000 ) / 2;                        % causality imposed
end
% (1) compute the hollowed median filter predictor (Eq. 3.3)
[ ~, pred ]                     = cch_conv( cc( :, u1, u2 ), W, convType, hollowF, 0 );

% (2) compute conditional rate CCH (Eq. 4)
den                             =  nspks1  * dt; 
crCCH                           = ( cc( :, u1, u2 ) - pred ) ./ den;      	% [spk/s]

% (3) compute the STGs (Eq. 5)
t                               = cchbins * 1000;                           % [ms]
t_ROI                           = t >= roiMS( 1 ) & t <= roiMS( 2 );
g1                              = calc_stg( crCCH, t_ROI, dt );

% (4) plot CCH and dcCCH for the selected pair
subplot( 2, 2, 3 );
plot_CCH( cchbins, cc( :, u1, u2 ), Cbars, NaN, g1 )

subplot( 2, 2, 4 );
cidx                            = ( u2 - 1 ) * bs + u1;
plot_CCH( cchbins, dcCCH( :, cidx ), Cdc, NaN, eSTG1( cidx ) )

return

%------------------------------------------------------------------------
% LOCAL FUNCTIONS
function plot_CCH( cchbins, cch, ColorBar, rSTG, eSTG )                     % plots CCH as a bar graph

bar( cchbins, cch, 1, 'FaceColor', ColorBar, 'EdgeColor', 'none' );
set( gca, 'box', 'off', 'tickdir', 'out' )
xlabel( 'Time lag [s]' )
ylabel( 'Count' )
title( sprintf( 'rSTG: %0.5f  eSTG: %0.5f', [ rSTG eSTG ] ) )

return

% EOF

