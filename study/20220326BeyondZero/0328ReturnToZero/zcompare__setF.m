if ( -1 == fSeed )
	fSeed = setprngstates();
else
	setprngstates(fSeed);
endif
%
switch( fType )
case 10
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = eye(sizeF,sizeX) + 1.0e-1*randn(sizeF,sizeX);
	matA0 = 1.0e-3*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-3*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 50
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = eye(sizeF,sizeX) + randn(sizeF,sizeX);
	matA0 = 1.0e-2*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-2*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
otherwise
	error( "Invalid fType." );
endswitch
%
runFStr = sprintf( "%d(%d)x%d", fType, fSeed, sizeX );
msg( __FILE__, __LINE__, sprintf( "Generated F '%s'.", runFStr ) );
