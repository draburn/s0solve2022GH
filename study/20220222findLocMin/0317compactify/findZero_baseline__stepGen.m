function [ vecX_next, vecF_next, stepGen_datOut ] = findZero_baseline__stepGen( vecX, vecF, funchXTrialOfP, funchFModel, funchF, stepGen_prm )
	stepGen_datOut = [];
	stepGen_datOut.fevalCount = 0;
	%
	switch (tolower(mygetfield( stepGen_prm, "stepType", "geombt"  )))
	case { "blind", "full" }
		p = 1.0;
		vecX_next = funchXTrialOfP(p);
		vecF_next = funchF(vecX_next); stepGen_datOut.fevalCount++;
	case { "geombt" }
		p = 1.0;
		fTol = mygetfield( stepGen_prm, "fTol", eps );
		xTol = mygetfield( stepGen_prm, "xTol", eps );
		btFactor = mygetfield( stepGen_prm, "btFactor", 3.0 );
		for emergencyBreak=1:1E8
			assert( emergencyBreak < 1000 );
			vecX_next = funchXTrialOfP(p);
			vecF_next = funchF(vecX_next); stepGen_datOut.fevalCount++;
			if ( norm(vecF_next) < norm(vecF) - fTol )
				break;
			elseif ( norm(vecX_next-vecX) < xTol )
				msg( __FILE__, __LINE__, "Reached backtracking limit." );
				break;
			endif
			p /= btFactor;
		endfor
	case { "minscan" }
		funchOmegaOfP = @(p)( sumsq(funchF(funchXTrialOfP(p)))/2.0 );
		fminbnd_prm = optimset( "TolX", eps^2, "TolFun", eps^2 );
		[ p, fminbnd_fval, fminbnd_info, fminbnd_output ] = fminbnd( funchOmegaOfP, 0.0, 1.0, fminbnd_prm );
		stepGen_datOut.fevalCount += fminbnd_output.funcCount;
		vecX_next = funchXTrialOfP(p);
		vecF_next = funchF(vecX_next); stepGen_datOut.fevalCount++;
	otherwise
		error( "Invalid value of stepType." );
	endswitch
	stepGen_datOut.p = p;
return;
endfunction
