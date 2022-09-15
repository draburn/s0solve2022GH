function xDatOut = groot_x( funchF, vecX0, solveSetPrm=[], default_solverPrm=[], xPrm=[] )
	startTime = time();
	if ( isempty(solveSetPrm) )
		solveSetPrm.n = 2;
		solveSetPrm.s(1).f = @groot_fsolve;
		solveSetPrm.s(2).f = @groot_jfnk_basic;
	endif
	%
	sizeX = size( vecX0, 1 );
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = funchF( vecX0 );
	sizeF = size( vecF0, 1 );
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	assert( isposintscalar(solveSetPrm.n) );
	assert( isvector(solveSetPrm.s) );
	assert( max(size(solveSetPrm.s)) == solveSetPrm.n );
	for n = 1 : solveSetPrm.n
		foo.f = solveSetPrm.s(n).f;
		assert( is_function_handle( foo.f ) );
		foo.fInfo = functions( foo.f );
		foo.strSolverName = trimFileName(foo.fInfo.function)(7:end);
		%
		foo.solverPrm = mygetfield( solveSetPrm.s(n), "p", [] );
		foo.solverPrm = overwritefields( default_solverPrm, foo.solverPrm );
		%
		foo.startTime = time();
		[ foo.vecXBest, foo.strGrootFlag, foo.fevalCount, foo.datOut ] = foo.f( funchF, vecX0, foo.solverPrm );
		foo.elapsedTime = time()-foo.startTime();
		assert( isrealarray(foo.vecXBest,[sizeX,1]) );
		assert( ischar(foo.strGrootFlag) );
		assert( isvector(foo.strGrootFlag) );
		assert( isposintscalar(foo.fevalCount) );
		%
		foo.vecFBest = funchF( foo.vecXBest );
		assert( isrealarray(foo.vecFBest,[sizeF,1]) );
		foo.fBest = norm(foo.vecFBest);
		%
		msg( __FILE__, __LINE__, sprintf( "  %10s:  %7s  %7d  %10.3e", foo.strSolverName, foo.strGrootFlag, foo.fevalCount, foo.fBest ) );
		%
		xDatOut.s(n).strGrootFlag = foo.strGrootFlag;
		xDatOut.s(n).fevalCount = foo.fevalCount;
		xDatOut.s(n).fBest = foo.fBest;
		xDatOut.s(n).matInfoA = foo.datOut.matInfoA;
		xDatOut.s(n).matInfoB = foo.datOut.matInfoB;
		xDatOut.s(n).elapsedTime = foo.elapsedTime;
		%
		clear foo;
	endfor
	xDatOut.elapsedTime = time() - startTime;
return;
endfunction
