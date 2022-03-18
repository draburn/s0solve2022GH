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
	otherwise
		error( "Invalid value of stepType." );
	endswitch
	stepGen_datOut.p = p;
return;
endfunction
