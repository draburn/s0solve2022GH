% Function...
%  A 'gnostic for studying how various forms of Jacobian updating impact convergence.

function [ vecX, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchF, prm=[] )
	commondefs;
	findLocMin_gnostic_jupdate__defs;
	time0 = time();
	fevalCount = 0;
	jevalCount = 0;
	collected_vecXVals = [];
	collected_vecFVals = [];
	%
	verbLev = mygetfield( prm, "verbLev", VERBLEV__COPIOUS );
	valdLev = mygetfield( prm, "valdLev", VALDLEV__UNLIMITED );
	if ( valdLev >= VALLEV__LOW )
		assert( isrealscalar(verbLev) );
		assert( isrealscalar(valdLev) );
	endif
	%
	sizeX = size(vecX0,1);
	if ( valdLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecX0,[sizeX,1]) );
	endif
	%
	vecF0 = funchF( vecX0 );
	fevalCount++;
	sizeF = size(vecF0,1);
	collected_vecXVals = [ collected_vecXVals, vecX0 ];
	collected_vecFVals = [ collected_vecFVals, vecF0 ];
	if ( valdLev >= VALLEV__MEDIUM )
		assert( isrealarray(vecF0,[sizeF,1]) );
	endif
	omega0 = sumsq(vecF0)/2.0;
	%
	matJ0 = jacobs( vecX0, funchF ); jevalCount++;
	if ( valdLev >= VALLEV__MEDIUM )
		assert( isrealarray(matJ0,[sizeF,sizeX]) );
	endif
	%
	%
	%
	iterMax = mygetfield( prm, "iterMax", 300 );
	if ( valdLev >= VALLEV__LOW )
		assert( isrealscalar(iterMax) );
	endif
	stepType = mygetfield( prm, "stepType", STEP_TYPE__BLIND_NEWTON );
	if ( valdLev >= VALLEV__LOW )
		assert( isrealscalar(stepType) );
	endif
	%
	jupdateType = mygetfield( prm, "jupdateType", JUPDATE_TYPE__RECALC );
	if ( valdLev >= VALLEV__LOW )
		assert( isrealscalar(jupdateType) );
	endif
	%
	cholSafeTol = mygetfield( prm, "cholSafeTol", sqrt(eps) );
	if ( valdLev >= VALLEV__LOW )
		assert( isrealscalar(cholSafeTol) );
		assert( 0.0 < cholSafeTol );
	endif
	%
	%
	%	
	if ( nargout >= 2 )
		datOut = [];
	endif
	%
	%
	%
	vecX = vecX0;
	vecF = vecF0;
	matJ = matJ0;
	iterCount = 0;
	while (1)
		omega = sumsq(vecF)/2.0;
		%
		if ( verbLev >= VERBLEV__PROGRESS )
		if ( abs( iterCount - round(sqrt(iterCount))^2 ) < 0.001 )
		if ( 0 == iterCount )
			msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d,  %3d;  %10.3e, %10.3e.", ...
			  time()-time0, iterCount, ...
			  fevalCount, jevalCount, ...
			  sumsq(vecF)/2.0, -1.0 ) );
		else
			msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d,  %3d;  %10.3e, %10.3e.", ...
			  time()-time0, iterCount, ...
			  fevalCount, jevalCount, ...
			  sumsq(vecF)/2.0, (sumsq(vecF_prev)-sumsq(vecF))/2.0 ) );
		endif
		endif
		endif
		%
		if ( 0.0 == norm(vecF) )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "Reached vecF = 0.0." );
			break;
		elseif ( omega < eps^2*omega0 )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "Reached omega < eps^2*omega0." );
			break;
		elseif ( 0 < iterCount )
			if ( reldiff( vecF, vecF_prev ) == 0.0 )
				msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "There was no change in vecF." );
				break;
			endif
		endif
		%
		iterCount++;
		if ( iterCount > iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "Reached iterMax." );
			break;
		endif
		%
		if ( valdLev >= VALDLEV__MEDIUM )
			assert( isrealarray(vecX,[sizeX,1]) );
			assert( isrealarray(vecF,[sizeF,1]) );
			assert( isrealarray(matJ,[sizeF,sizeX]) );
		endif
		%
		%
		%
		switch( stepType )
		case STEP_TYPE__BLIND_NEWTON
			vecG = matJ'*vecF;
			matJTJ = matJ'*matJ;
			[ matR, cholFlag ] = chol( matJTJ );
			if ( 0~=cholFlag || min(diag(matR)) <= cholSafeTol*max(abs(diag(matR))) )
				hScale = max(abs(diag(matJTJ)));
				assert( 0.0 ~= hScale );
				[ matR, cholFlag ] = chol( matJ'*matJ + sqrt(eps)*hScale*eye(sizeX,sizeX) );
				if ( 0~=cholFlag || min(diag(matR)) <= cholSafeTol*max(abs(diag(matR))) )
					error( "Bad!" );
				endif
			endif
			msgif( verbLev >= VERBLEV__UNLIMITED, __FILE__, __LINE__, "Evaluating vecDelta = matR \ ( matR' \ (-vecG) )." );
			vecDelta = matR \ ( matR' \ (-vecG) );
			if ( valdLev >= VALDLEV__MEDIUM )
				assert( isrealarray(vecDelta,[sizeX,1]) );
			endif
			msgif( verbLev >= VERBLEV__UNLIMITED, __FILE__, __LINE__, "Evaluated vecDelta = matR \ ( matR' \ (-vecG) )." );
			vecX_next = vecX + vecDelta;
			vecF_next = funchF( vecX_next );
			fevalCount++;
			if ( valdLev >= VALDLEV__MEDIUM )
				assert( isrealarray(vecF_next,[sizeF,1]) );
			endif
		case STEP_TYPE__BLIND_GRAD_MIN
			error( "STEP_TYPE__BLIND_GRAD_MIN is not implemented yet." );
		case STEP_TYPE__SCAN_LEV_MIN
			error( "STEP_TYPE__SCAN_LEV_MIN is not implemented yet." );
		otherwise
			error( "Invalid value of stepType." );
		endswitch
		%
		%
		%
		switch( jupdateType )
		case JUPDATE_TYPE__NONE
			matJ_next = matJ;
		case JUPDATE_TYPE__BROYDEN
			fooX = vecX_next - vecX;
			fooY = vecF_next - ( vecF + matJ*fooX );
			if ( 0.0 == fooX'*fooX )
				msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "Step size was zero." );
				break;
			endif
			fooJ = (fooY*(fooX'))/(fooX'*fooX);
			matJ_next = matJ + fooJ;
			if ( valdLev >= VALLEV__HIGH )
				assert( reldiff( matJ_next*(vecX_next-vecX), vecF_next-vecF ) <= sqrt(eps) );
			endif
			clear fooX;
			clear fooY;
			clear fooJ;
		case JUPDATE_TYPE__SECANT_REORTHONORM
			error( "JUPDATE_TYPE__SECANT_REORTHONORM is not implemented yet." );
		case JUPDATE_TYPE__LESQUJ_PRIMAL
			collected_vecXVals = [ collected_vecXVals, vecX_next ];
			collected_vecFVals = [ collected_vecFVals, vecF_next ];
			lesquj_prm = [];
			lesquj_prm.jevalDat(1).vecX = vecX0;
			lesquj_prm.jevalDat(1).vecF = vecF0;
			lesquj_prm.jevalDat(1).matJ = matJ0;
			%lesquj_prm.useDistanceWeights = false; % Because caues error.
			leqsuj_prm.useLatestPtAs0 = true;
			msgif( verbLev >= VERBLEV__UNLIMITED, __FILE__, __LINE__, "Calling calcLesquj_basic()." );
			[ lesquj_vecX0, lesquj_vecF0, lesquj_matJ0, lesqu_datOut ] = calcLesquj_basic( collected_vecXVals, collected_vecFVals, lesquj_prm );
			msgif( verbLev >= VERBLEV__UNLIMITED, __FILE__, __LINE__, "Back from calcLesquj_basic()." );
			matJ_next = lesquj_matJ0;
			clear lesquj_vecX0;
			clear lesquj_vecF0;
			clear lesqu_datOut;
			clear lesquj_prm;
		case JUPDATE_TYPE__RECALC
			matJ_next = jacobs( vecX, funchF ); jevalCount++;
		otherwise
			error( "Invalid value of jupdateType." );
		endswitch
		%
		%
		%
		vecX_prev = vecX;
		vecF_prev = vecF;
		matJ_prev = matJ;
		%
		vecX = vecX_next;
		vecF = vecF_next;
		matJ = matJ_next;
	endwhile
	%
	if ( verbLev >= VERBLEV__MAIN )
	if ( 0 == iterCount )
		msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d,  %3d;  %10.3e, %10.3e.", ...
		  time()-time0, iterCount, ...
		  fevalCount, jevalCount, ...
		  sumsq(vecF)/2.0, -1.0 ) );
	else
		msg( __FILE__, __LINE__, sprintf( "  %10.3e, %4d;  %5d,  %3d;  %10.3e, %10.3e.", ...
		  time()-time0, iterCount, ...
		  fevalCount, jevalCount, ...
		  sumsq(vecF)/2.0, (sumsq(vecF_prev)-sumsq(vecF))/2.0 ) );
	endif
	endif
endfunction

%!function [ vecF, matJ ] = funcFJ_cubyDiagTest( vecX, c )
%!	sizeX = size(vecX,1);
%!	vecXE = (1:sizeX)';
%!	vecF = (vecX-vecXE) + c*(vecX-vecXE).^3;
%!	if ( nargout >= 2 )
%!		matJ = eye(sizeX) + c*3.0*diag((vecX-vecXE).^2);
%!	endif
%!endfunction


%!test
%!	clear;
%!	commondefs;
%!	findLocMin_gnostic_jupdate__defs;
%!	setprngstates(0);
%!	numFigs = 0;
%!	%
%!	caseNum = 40;
%!	msg( __FILE__, __LINE__, sprintf( "caseNum = %d.", caseNum ) );
%!	switch (caseNum)
%!	case 0
%!		sizeX = 2;
%!		c_cuby = 0.0;
%!	case 10
%!		sizeX = 15;
%!		c_cuby = 0.0;
%!	case 20
%!		sizeX = 2;
%!		c_cuby = 0.01;
%!	case 30
%!		sizeX = 15;
%!		c_cuby = 0.01;
%!	case 40
%!		sizeX = 15;
%!		c_cuby = 0.1;
%!	case 100
%!		sizeX = 15;
%!		c_cuby = 1.0;
%!	case 200
%!		sizeX = 20;
%!		c_cuby = 1.0;
%!	otherwise
%!		error( "Ivalid caseNum." );
%!	endswitch
%!	funchFJ = @(dummyX)( funcFJ_cubyDiagTest( dummyX, c_cuby ) );
%!	vecX0 = zeros(sizeX,1);
%!	%
%!	msg( __FILE__, __LINE__, "" );
%!	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__NONE ~~~ " );
%!	prm = [];
%!	prm.jupdateType = JUPDATE_TYPE__NONE;
%!	[ vecXF, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchFJ, prm );
%!	vecX0 = zeros(sizeX,1);
%!	%
%!	msg( __FILE__, __LINE__, "" );
%!	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__BROYDEN ~~~ " );
%!	prm = [];
%!	prm.jupdateType = JUPDATE_TYPE__BROYDEN;
%!	[ vecXF, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchFJ, prm );
%!	%
%!	msg( __FILE__, __LINE__, "" );
%!	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__LESQUJ_PRIMAL ~~~ " );
%!	prm = [];
%!	prm.jupdateType = JUPDATE_TYPE__LESQUJ_PRIMAL;
%!	[ vecXF, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchFJ, prm );
%!	vecX0 = zeros(sizeX,1);
%!	%
%!	msg( __FILE__, __LINE__, "" );
%!	msg( __FILE__, __LINE__, "~~~ JUPDATE_TYPE__RECALC ~~~ " );
%!	prm = [];
%!	prm.jupdateType = JUPDATE_TYPE__RECALC;
%!	[ vecXF, datOut ] = findLocMin_gnostic_jupdate( vecX0, funchFJ, prm );
