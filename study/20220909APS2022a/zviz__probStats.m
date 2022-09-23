	bigNVals = zeros(1,numProbs);
	condVals = zeros(1,numProbs);
	for n = 1 : numProbs
		this_probType = zcdo.prob(n).probType;
		this_sizeXIsh = zcdo.prob(n).sizeXIsh;
		this_probSeed = zcdo.prob(n).probSeed;
		this_probGenPrm = zcdo.prob(n).probGenPrm;
		[ this_funchF, this_vecX0, this_probDat ] = genFunchAPS2022_fromType( this_probType, this_sizeXIsh, this_probSeed, this_probGenPrm );
		%
		foo = this_probDat.genFunchDatOut;
		this_matJ0Secret = foo.matA11L*foo.matA11R + foo.matA12L*foo.matA12R + foo.matB1;
		%
		this_probStr = sprintf( "Prob %d", n );
		this_probSizeStr = sprintf( "(%dx%d)", this_probDat.sizeF, this_probDat.sizeX );
		this_cond = cond(this_matJ0Secret);
		%
		bigNVals(n) = this_probDat.sizeX;
		condVals(n) = this_cond;
		if (0)
			msg( __FILE__, __LINE__, sprintf( ...
			  "  %12s:  %12s, cond0 %10.3e", ...
			  this_probStr, this_probSizeStr, this_cond ) );
		endif
	endfor
	clear n;
	clear this*;
	clear foo;
	%
	msg( __FILE__, __LINE__, sprintf( ...
	  "Size range: %d ~ %d; cond range: %g ~ %g.", ...
	  min(bigNVals), max(bigNVals), min(condVals), max(condVals) ) );
