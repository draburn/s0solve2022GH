function [ vecX_step, vecF_step, datOut ] = findZero_baseline__step( vecX, funchDeltaOfP, funchF, prm );
	pLo = mygetfield( prm, "pLo", 0.0 );
	pHi = mygetfield( prm, "pHi", 1.0 );
	vecDelta = funchDeltaOfP(pHi);
	vecX_step = vecX + vecDelta;
	vecF_step = funchF( vecX_step );
	datOut = [];
	return
	error( "CODE BELOW IS JUNK." );
		%
		%
		stepSizeMode = mygetfield( prm, "stepSizeMode", "tr" );
		switch ( tolower(stepSizeMode) )
		case { "blind", "full" }
			p = 1.0;
			vecDelta = funchDeltaOfP(p);
			vecX_trial = vecX + vecDelta;
			vecDelta = vecX_trial - vecX;
			%vecF_trial = 
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
