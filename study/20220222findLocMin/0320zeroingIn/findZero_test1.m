	clear;
	setVerbLevs;
	setprngstates(14507280); % 14507280: pool happens to work better.
	numFigs = 0;
	%
	sizeX = 30;
	sizeF = 30;
	vecXE = randn(sizeX,1);
	matJE = randn(sizeF,sizeX);
	matA0 = 0.001*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 0.001*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	%
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = zeros(sizeX,1);
	%
	[ vecXF_fg, vecFF_fg, datOut_fg ] = findZero_fsolveGnostic( vecX0, funchF );
	[ vecXF_basic, vecFF_basic, datOut_basic ] = findZero_basic( vecX0, funchF );
	[ vecXF_mimic, vecFF_mimic, datOut_mimic ] = findZero_mimic( vecX0, funchF );
	[ vecXF_beyond, vecFF_beyond, datOut_beyond ] = findZero_beyond( vecX0, funchF );
	%
	numFigs++; figure( numFigs );
	semilogy( ...
	  datOut_fg.fevalCountVals, datOut_fg.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
	  datOut_basic.fevalCountVals, datOut_basic.fNormVals+eps, 'x-', 'markersize', 20, 'linewidth', 2, ...
	  datOut_mimic.fevalCountVals, datOut_mimic.fNormVals+eps, 's-', 'markersize', 20, 'linewidth', 2, ...
	  datOut_beyond.fevalCountVals, datOut_beyond.fNormVals+eps, 'p-', 'markersize', 20, 'linewidth', 2 );
	grid on;
	ylabel( "||f||" );
	xlabel( "feval count" );
	%
	numFigs++; figure( numFigs );
	semilogy( ...
	  datOut_fg.iterCountVals, datOut_fg.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
	  datOut_basic.iterCountVals, datOut_basic.fNormVals+eps, 'x-', 'markersize', 20, 'linewidth', 2, ...
	  datOut_mimic.iterCountVals, datOut_mimic.fNormVals+eps, 's-', 'markersize', 20, 'linewidth', 2, ...
	  datOut_beyond.iterCountVals, datOut_beyond.fNormVals+eps, 'p-', 'markersize', 20, 'linewidth', 2 );
	grid on;
	ylabel( "||f||" );
	xlabel( "iteration count" );
