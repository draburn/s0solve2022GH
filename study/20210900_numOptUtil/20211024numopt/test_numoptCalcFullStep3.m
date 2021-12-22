	clear;
	commondefs;
	thisFile = "test_numoptCalcFullStep1";
	msg( thisFile, __LINE__, "Starting test." );
	setprngstates(0);
	tic();
	%
	msg( thisFile, __LINE__, "This is a positive definite case without a root." );
	sizeX = 2;
	sizeF = 3;
	vecFFoo = randn(sizeF,1);
	matJ = randn(sizeF,sizeX);
	funchF = @(x)( vecFFoo + matJ*x );
	vecX0 = zeros(sizeX,1);
	vecF0 = funchF( vecX0 );
	%
	omega0 = 0.5 * vecF0' * vecF0;
	msg( thisFile, __LINE__, sprintf( "omega0 = %e.", omega0 ) );
	vecG = matJ' * vecF0;
	matH = matJ' * matJ;
	%
	vecDelta = numoptCalcFullStep( omega0, vecG, matH );
	vecX1 = vecX0 + vecDelta;
	vecF1 = funchF( vecX1 );
	omega1 = 0.5 * vecF1' * vecF1;
	msg( thisFile, __LINE__, sprintf( "omega1 = %e.", omega1 ) );
	assert( omega1 <= omega0 )
	%
	msg( thisFile, __LINE__, "Finished test." );
	toc();
