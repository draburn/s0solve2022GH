if ( -1 == fSeed )
	fSeed = setprngstates();
else
	setprngstates(fSeed);
endif
%
switch( fType )
case 1
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = eye(sizeF,sizeX) + 1.0e-2*randn(sizeF,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) );
	%
	vecX0 = randn(sizeX,1);
case 3
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = eye(sizeF,sizeX) + 1.0e-1*randn(sizeF,sizeX);
	matA0 = 1.0e-6*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-6*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	if (0)
	function f = foo(x,vecXE,matJE,matA0,matA1,matA2,matB0,matB1,matB2,matB3)
		msg( __FILE__, __LINE__, "feval!" );
		y = @(x)( x - vecXE );
		f = matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) );
	endfunction
	funchF = @(x)(foo(x,vecXE,matJE,matA0,matA1,matA2,matB0,matB1,matB2,matB3));
	endif
	%
	vecX0 = randn(sizeX,1);
case 5
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = eye(sizeF,sizeX) + 1.0e-2*randn(sizeF,sizeX);
	matA0 = 1.0e-4*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-4*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 6
	sizeF = 2*sizeX;
	vecXE = randn(sizeX,1);
	matJE = eye(sizeF,sizeX) + 1.0e-2*randn(sizeF,sizeX);
	matA0 = 1.0e-4*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-4*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 7
	sizeF = sizeX+1;
	vecXE = randn(sizeX,1);
	matJE = eye(sizeF,sizeX) + 1.0e-2*randn(sizeF,sizeX);
	matA0 = 1.0e-4*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-4*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
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
case 60
	sizeF = sizeX;
	vecXE = ones(sizeX,1);
	vecFE = ones(sizeF,1);
	funchF = @(x)( vecFE + (x-vecXE).^2 );
	vecX0 = zeros(sizeX,1);
case 90
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = randn(sizeF,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) );
	%
	vecX0 = randn(sizeX,1);
	matJA = matJE;
case 91
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = randn(sizeF,sizeX);
	matA0 = 1.0e-8*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-8*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
	matJA = (1.0+1.0e-1*randn(sizeF,sizeX)).*matJE + 1.0e-1*randn(sizeF,sizeX);
case 300
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = ones(sizeF,sizeX);
	matA0 = 0.0e0*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 0.0e0*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 313
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = diag(1.0+0.1*abs(randn(sizeX,1))); + 1.0e-8*randn(sizeF,sizeX);
	matA0 = 1.0e-8*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-8*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 314
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = diag(1.0+0.1*abs(randn(sizeX,1))); + 1.0e-8*randn(sizeF,sizeX);
	matA0 = 1.0e-8*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-8*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	%
	function vecF = funcF( vecX, vecXE, matJE, matA0, matA1, matA2, matB0, matB1, matB2, matB3 )
		vecY = vecX - vecXE;
		vecF =  matJE*vecY + matA0*( (matA1*vecY) .* (matA2*vecY) ) + matB0*( (matB1*vecY) .* (matB2*vecY) .* (matB3*vecY) );
		omega = sumsq(vecF)/2.0;
		msg( __FILE__, __LINE__, sprintf( "funchF: omega = %0.6e.", omega ) );
		return;
	endfunction
	funchF = @(x)( funcF( x, vecXE, matJE, matA0, matA1, matA2, matB0, matB1, matB2, matB3 ) );
	%
	vecX0 = randn(sizeX,1);
	feval_vecXVals(:,1) = vecX0;
	%save( "feval_vecXVals.m", "feval_vecXVals" );
case 319
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = diag(1:sizeX); + 1.0e-8*randn(sizeF,sizeX);
	matA0 = 1.0e-8*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-8*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 320
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = diag(1.0+abs(10.0*randn(sizeX,1))); + 1.0e-8*randn(sizeF,sizeX);
	matA0 = 1.0e-8*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-8*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 510
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = diag(randn(sizeX)) + 1.0E-2*randn(sizeF,sizeX);
	matA0 = 1.0e-7*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-7*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 520
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = diag(randn(sizeX)) + 1.0E-1*randn(sizeF,sizeX);
	matA0 = 1.0e-6*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-6*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 530
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = diag(randn(sizeX)) + 1.0E-1*randn(sizeF,sizeX);
	matA0 = 1.0e-4*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-4*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 540
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = randn(sizeF,sizeX);
	matA0 = randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 1000
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = randn(sizeF,sizeX);
	matA0 = 0.0e0*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 0.0e0*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 1010
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = randn(sizeF,sizeX);
	matA0 = 1.0e-6*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-6*randn(sizeF,sizeX);
	matB1 = randn(sizeX,sizeX);
	matB2 = randn(sizeX,sizeX);
	matB3 = randn(sizeX,sizeX);
	y = @(x)( x - vecXE );
	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
	%
	vecX0 = randn(sizeX,1);
case 1020
	sizeF = sizeX;
	vecXE = randn(sizeX,1);
	matJE = randn(sizeF,sizeX);
	matA0 = 1.0e-5*randn(sizeF,sizeX);
	matA1 = randn(sizeX,sizeX);
	matA2 = randn(sizeX,sizeX);
	matB0 = 1.0e-5*randn(sizeF,sizeX);
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
