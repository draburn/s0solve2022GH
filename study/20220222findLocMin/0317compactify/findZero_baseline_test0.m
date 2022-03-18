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
	[ vecXF, vecFF, datOut ] = findZero_baseline( vecX0, funchF )
