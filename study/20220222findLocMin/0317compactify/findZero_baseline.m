function [ vecXF, vecFF, datOut ] = findZero_baseline( vecX0, funchF, prm=[] )
	time0 = time();
	fevalCount = 0;
	setVerbLevs;
	verbLev = mygetfield( prm, "verbLev", VERBLEV__UNLIMITED );
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = funchF( vecX0 ); fevalCount++;
	sizeF = size(vecF0,1);
	assert( isrealarray(vecF0,[sizeF,1]) );
	%
	%
	bestIsNot0 = false;
	vecX_best = vecX0;
	vecF_best = vecF0;
	datOut = [];
	%
	%
	vecX_prev = [];
	vecF_prev = [];
	matJ_prev = [];
	vecX = vecX0;
	vecF = vecF0;
	iterCount = 0;
	while (1)
		iterCount++;
		fTol = mygetfield( prm, "fTol", eps );
		iterMax = mygetfield( prm, "iterMax", 50 );
		if ( norm(vecF) <= fTol )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "STRONG SUCCESS: norm(vecF) <= fTol." );
			break;
		elseif ( iterCount >= iterMax )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: iterCount >= iterMax." );
			break;
		endif
		%
		%
		%
		modelGen_prm = mygetfield( prm, "modelGen_prm", [] );
		[ matJ, ary3Kappa, funchFModel, modelGen_datOut ] = findZero_baseline__modelGen( vecX, vecF, vecX_prev, vecF_prev, matJ_prev, funchF, modelGen_prm );
		fevalCount += modelGen_datOut.fevalCount;
		%
		%
		%
		curveGen_prm = mygetfield( prm, "curveGen_prm", [] );
		[ funchXTrialOfP, curveGen_datOut ] = findZero_baseline__curveGen( vecX, vecF, matJ, ary3Kappa, funchF, curveGen_prm );
		fevalCount += curveGen_datOut.fevalCount; % Always zero?
		%
		%
		%
		stepGen_prm = mygetfield( prm, "stepGen_prm", [] );
		[ vecX_next, vecF_next, stepGen_datOut ] = findZero_baseline__stepGen( vecX, vecF, funchXTrialOfP, funchFModel, funchF, stepGen_prm );
		fevalCount += stepGen_datOut.fevalCount;
		%
		%
		%
		if ( norm(vecF_next) < norm(vecF_best) )
			bestIsNot0 = true;
			vecX_best = vecX_next;
			vecF_best = vecF_next;
		endif
		if ( norm(vecF_next) >= norm(vecF) )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "ALGORITHM BREAKDOWN: norm(vecF_next) >= norm(vecF)." );
			break;
		elseif ( norm(vecF_next) >= norm(vecF) - fTol )
			msgif( verbLev >= VERBLEV__MAIN, __FILE__, __LINE__, "IMPOSED STOP: norm(vecF_next) >= norm(vecF) - fTol." );
			break;
		endif
		%
		vecX_prev = vecX;
		vecF_prev = vecF;
		matJ_prev = matJ;
		vecX = vecX_next;
		vecF = vecF_next;
		%
		msgif( verbLev >= VERBLEV__PROGRESS, __FILE__, __LINE__, sprintf( "  %10.3e, %4d, %5d;  %10.3e;  %10.3e, %10.3e, %10.3e;  %10.3e, %10.3e, %10.3e.", ...
		  time()-time0, iterCount, fevalCount, ...
		  norm(matJ_prev'*vecF_prev), ...
		  norm(vecX-vecX0), norm(vecX-vecX0)-norm(vecX_prev-vecX0), norm(vecX-vecX_prev), ...
		  norm(vecF), norm(vecF_prev)-norm(vecF), norm(vecF_prev-vecF) ) );	
	endwhile
	vecXF = vecX_best;
	vecFF = vecF_best;
	datOut.fevalCount = fevalCount;
return;
endfunction


%!test
%!	setprngstates(0);
%!	%
%!	sizeX = 3;
%!	sizeF = 3;
%!	%
%!	vecXE = randn(sizeX,1);
%!	matJE = randn(sizeF,sizeX);
%!	matA0 = 0.1*randn(sizeF,sizeX);
%!	matA1 = randn(sizeX,sizeX);
%!	matA2 = randn(sizeX,sizeX);
%!	matB0 = 0.1*randn(sizeF,sizeX);
%!	matB1 = randn(sizeX,sizeX);
%!	matB2 = randn(sizeX,sizeX);
%!	matB3 = randn(sizeX,sizeX);
%!	y = @(x)( x - vecXE );
%!	funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
%!	%
%!	vecX0 = zeros(sizeX,1);
%!	%
%!	[ vecXF, vecFF, datOut ] = findZero_baseline( vecX0, funchF )