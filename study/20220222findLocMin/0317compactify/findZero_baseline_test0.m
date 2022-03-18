	clear;
	setprngstates(0);
	numFigs = 0;
	%
	sizeX = 3;
	sizeF = 3;
	vecXE = randn(sizeX,1);
	matJE = randn(sizeF,sizeX);
	matA0 = 0.1*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 0.1*randn(sizeF,sizeX);
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
	%
	prm_broyd = [];
	prm_broyd.modelGen_prm.useInexactJ = "true";
	[ vecXF_broyd, vecFF_broyd, datOut_broyd ] = findZero_baseline( vecX0, funchF, prm_broyd );
	%
	[ vecXF, vecFF, datOut ] = findZero_baseline( vecX0, funchF );
	%
	numFigs++; figure( numFigs );
	semilogy( ...
	  datOut_fg.fevalCountVals, datOut_fg.fNormVals+eps, 'o-', 'markersize', 20, 'linewidth', 2, ...
	  datOut_broyd.fevalCountVals, datOut_broyd.fNormVals+eps, 'x-', 'markersize', 20, 'linewidth', 2, ...
	  datOut.fevalCountVals, datOut.fNormVals+eps, 'p-', 'markersize', 20, 'linewidth', 2 );
	grid on;
