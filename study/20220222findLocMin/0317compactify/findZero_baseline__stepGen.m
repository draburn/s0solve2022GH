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
		fTol = mygetfield( stepGen_prm, "fTol", eps^0.7 );
		xTol = mygetfield( stepGen_prm, "xTol", eps^0.7 );
		emergencyBreak = 0;
		while (1)
			vecX_next = funchXTrialOfP(p);
			vecF_next = funchF(vecX_next); stepGen_datOut.fevalCount++;
			if ( norm(vecF_next) < norm(vecF) - fTol )
				break;
			elseif ( norm(vecX_next-vecX) < xTol )
				msg( __FILE__, __LINE__, "Reached backtracking limit." );
				break;
			endif
			p /= 3.0;
			emergencyBreak++;
			assert( emergencyBreak < 10000 );
		endwhile
	otherwise
		error( "Invalid value of stepType." );
	endswitch
	stepGen_datOut.p = p;
return;
endfunction
