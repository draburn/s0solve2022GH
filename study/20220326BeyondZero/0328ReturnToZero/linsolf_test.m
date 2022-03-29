	clear;
	setprngstates(0);
	numFigs = 0;
	%
	sizeX = 10;
	sizeF = 10;
	matA = eye(sizeF,sizeX) + 0.1 * randn(sizeF,sizeX);
	rcond(matA)
	vecX0 = zeros(sizeX,1);
	vecB = randn(sizeF,1);
	funchMatAProd = @(v)( matA*v );
	[ vecX, datOut ] = linsolf( funchMatAProd, vecB, vecX0 );
	size(datOut.matV,2)
	norm( (matA*vecX)-vecB ) / norm(vecB)
