function [ vecX_step, vecF_step, datOut ] = findZero_baseline__step( vecX, funchDeltaOfP, funchF, prm );
	doCurveViz = false;
	if (doCurveViz)
		findZero_baseline__curveViz;
		error( "Finizhed viz." );
	endif
	%
	pLo = mygetfield( prm, "pLo", 0.0 );
	pHi = mygetfield( prm, "pHi", 1.0 );
	%
	%
	datOut = [];
	datOut.fevalCount = 0;
	stepSizeMode = mygetfield( prm, "stepSizeMode", "full" );
	switch ( tolower(stepSizeMode) )
	case { "blind", "full" }
		p = pHi;
		vecDelta = funchDeltaOfP(p);
		vecX_step = vecX + vecDelta;
		vecF_step = funchF( vecX_step ); datOut.fevalCount++;
	case { "tr", "trust region" }
		error( "Not implemented." );
	case { "bt", "backtracking" }
		error( "Not implemented." );
	case { "scan" }
		funchOmegaOfP = @(p)( sumsq(funchF( vecX + funchDeltaOfP(p) ))/2.0 );
	otherwise
		error( "Invalid value of stepSizeMode." );
	endswitch
return;
endfunction
