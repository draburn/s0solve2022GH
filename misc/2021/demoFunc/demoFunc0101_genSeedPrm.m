function seedPrm = demoFunc0101_genSeedPrm( seedName = "trivial" )
	sizeX = 3;
	sizeF = sizeX;
	%
	coeffSeed_zero.c0 = 0.0;
	coeffSeed_zero.c1 = 0.0;
	coeffSeed_zero.c2 = 0.0;
	coeffSeed_fairVar.c0 = 0.0;
	coeffSeed_fairVar.c1 = 3.0;
	coeffSeed_fairVar.c2 = 3.0;
	%
	seedPrm.sizeX = sizeX;
	seedPrm.sizeF = sizeF;
	%
	%Use seedPrm.randState = floor(1e6*time) for a random state.
	seedPrm.randState = 0; % Should be applied to all relevant memebrs of rand() family.
	seedPrm.coeffSeed_x0     = coeffSeed_zero;
	seedPrm.coeffSeed_m1Diag = coeffSeed_zero;
	seedPrm.coeffSeed_m1Full = coeffSeed_zero;
	seedPrm.coeffSeed_m2Diag = coeffSeed_zero;
	seedPrm.coeffSeed_m2Full = coeffSeed_zero;
	seedPrm.coeffSeed_m3Diag = coeffSeed_zero;
	seedPrm.coeffSeed_m3Full = coeffSeed_zero;
	% Diagonal elements are sum of "full" and "diag" terms.
	% Much overwriting below.
	%
	if (strcmpi(seedName,"trivial"))
		seedPrm.coeffSeed_x0.c1 = 1.0;
		seedPrm.coeffSeed_m1Diag.c0 = 1.0;
		return;
	elseif (strcmpi(seedName,"diag-easy"))
		seedPrm.coeffSeed_x0.c1 = 1.0;
		seedPrm.coeffSeed_m1Diag.c0 = 1.0;
		seedPrm.coeffSeed_m1Diag.c1 = 0.5;
		return;
	elseif (strcmpi(seedName,"diag-moderate"))
		seedPrm.coeffSeed_x0.c1 = 5.0;
		seedPrm.coeffSeed_x0.c2 = 5.0;
		seedPrm.coeffSeed_m1Diag.c1 = 5.0;
		seedPrm.coeffSeed_m1Diag.c2 = 5.0;
		return;
	elseif (strcmpi(seedName,"lin-easy"))
		seedPrm.coeffSeed_x0.c1 = 1.0;
		seedPrm.coeffSeed_m1Diag.c0 = 1.0;
		seedPrm.coeffSeed_m1Diag.c1 = 0.3;
		seedPrm.coeffSeed_m1Full.c1 = 0.3 / sizeX;
		return;
	elseif (strcmpi(seedName,"lin-moderate"))
		seedPrm.coeffSeed_x0.c1 = 3.0;
		seedPrm.coeffSeed_x0.c2 = 3.0;
		seedPrm.coeffSeed_m1Diag.c0 = 1.0;
		seedPrm.coeffSeed_m1Diag.c1 = 0.4;
		seedPrm.coeffSeed_m1Full.c1 = 0.4 / sqrt(sizeX);
		return;
	elseif (strcmpi(seedName,"lin-tricky"))
		seedPrm.coeffSeed_x0.c1 = 5.0;
		seedPrm.coeffSeed_x0.c2 = 5.0;
		seedPrm.coeffSeed_m1Diag.c1 = 5.0;
		seedPrm.coeffSeed_m1Diag.c2 = 5.0;
		seedPrm.coeffSeed_m1Full.c1 = 1.0 / sqrt(sizeX);
		return;
	elseif (strcmpi(seedName,"lin-hard"))
		seedPrm.coeffSeed_x0.c1 = 5.0;
		seedPrm.coeffSeed_x0.c2 = 5.0;
		seedPrm.coeffSeed_m1Full.c1 = 5.0;
		seedPrm.coeffSeed_m1Full.c2 = 5.0;
		return;
	elseif (strcmpi(seedName,"easy"))
		seedPrm.coeffSeed_x0.c1 = 1.0;
		seedPrm.coeffSeed_m1Diag.c0 = 1.0;
		seedPrm.coeffSeed_m1Diag.c1 = 0.3;
		seedPrm.coeffSeed_m1Full.c1 = 0.3 / sizeX;
		seedPrm.coeffSeed_m2Diag.c0 = 0.3;
		seedPrm.coeffSeed_m2Diag.c1 = 0.1;
		seedPrm.coeffSeed_m2Full.c1 = 0.1 / sizeX;
		seedPrm.coeffSeed_m3Diag.c0 = 0.1;
		seedPrm.coeffSeed_m3Diag.c1 = 0.03;
		seedPrm.coeffSeed_m3Full.c1 = 0.03 / sizeX;
		return;
	elseif (strcmpi(seedName,"moderate"))
		seedPrm.coeffSeed_x0.c1 = 1.0;
		seedPrm.coeffSeed_m1Diag.c0 = 1.0;
		seedPrm.coeffSeed_m1Diag.c1 = 1.0;
		seedPrm.coeffSeed_m1Full.c1 = 0.05 / sizeX;
		%seedPrm.coeffSeed_m2Diag.c1 = 0.1;
		%seedPrm.coeffSeed_m2Full.c1 = 0.01 / sizeX;
		seedPrm.coeffSeed_m3Diag.c1 = 1.0;
		seedPrm.coeffSeed_m3Full.c1 = 0.01 / sizeX;
		return;
	elseif (strcmpi(seedName,"xhard"))
		seedPrm.coeffSeed_x0     = coeffSeed_fairVar;
		seedPrm.coeffSeed_m1Full = coeffSeed_fairVar;
		seedPrm.coeffSeed_m2Full = coeffSeed_fairVar;
		seedPrm.coeffSeed_m3Full = coeffSeed_fairVar;
		return;
	elseif (strcmpi(seedName,"bumpy"))
		seedPrm.coeffSeed_x0     = coeffSeed_fairVar;
		seedPrm.coeffSeed_m1Full = coeffSeed_fairVar;
		seedPrm.coeffSeed_m3Full = coeffSeed_fairVar;
		return;
	else
		error(sprintf( "Invalid seedName (\"%s\").\n", seedName ));
	endif
	%
return;
end

%!test
%!	seedPrm = demoFunc0101_genSeedPrm("easy");
