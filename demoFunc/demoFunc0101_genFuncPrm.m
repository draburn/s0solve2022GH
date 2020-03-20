function funcPrm = demoFunc0101_genFuncPrm( seedPrm )	
	randState_bkup  = rand("state");
	randnState_bkup  = randn("state");
	randeState_bkup  = rande("state");
	randgState_bkup  = randg("state");
	randpState_bkup  = randp("state");
	rand( "state", seedPrm.randState );
	%
	sizeX = seedPrm.sizeX;
	sizeF = seedPrm.sizeF;
	minSize = min([ sizeX, sizeF ]);
	%
	funcPrm.seedPrm = seedPrm;
	funcPrm.sizeX = sizeX;
	funcPrm.sizeF = sizeF;
	%
	% "Coefficient Generating Function".
	randSym = @(s1,s2)( (2.0*rand(s1,s2)) - 1.0 );
	cgf = @(s1,s2,c0,c1,c2)( c0 + (c1*randSym(s1,s2).*exp( c2 * randSym(s1,s2) )) );
	cgfWrap = @(s1,s2,cs)( cgf( s1, s2, cs.c0, cs.c1, cs.c2 ) );
	%
	funcPrm.x0  = cgfWrap( sizeX,   1,     seedPrm.coeffSeed_x0 );
	funcPrm.m1  = cgfWrap( sizeF,   sizeX, seedPrm.coeffSeed_m1Full ) ...
	      + diag( cgfWrap( minSize, 1,     seedPrm.coeffSeed_m1Diag ), sizeF, sizeX );
	funcPrm.m2a = cgfWrap( sizeF,   sizeX, seedPrm.coeffSeed_m2Full ) ...
	      + diag( cgfWrap( minSize, 1,     seedPrm.coeffSeed_m2Diag ), sizeF, sizeX );
	funcPrm.m2b = cgfWrap( sizeF,   sizeX, seedPrm.coeffSeed_m2Full ) ...
	      + diag( cgfWrap( minSize, 1,     seedPrm.coeffSeed_m2Diag ), sizeF, sizeX );
	funcPrm.m3a = cgfWrap( sizeF,   sizeX, seedPrm.coeffSeed_m3Full ) ...
	      + diag( cgfWrap( minSize, 1,     seedPrm.coeffSeed_m3Diag ), sizeF, sizeX );
	funcPrm.m3b = cgfWrap( sizeF,   sizeX, seedPrm.coeffSeed_m3Full ) ...
	      + diag( cgfWrap( minSize, 1,     seedPrm.coeffSeed_m3Diag ), sizeF, sizeX );
	funcPrm.m3c = cgfWrap( sizeF,   sizeX, seedPrm.coeffSeed_m3Full ) ...
	      + diag( cgfWrap( minSize, 1,     seedPrm.coeffSeed_m3Diag ), sizeF, sizeX );
	%
	rand("state",randState_bkup);
	randn("state",randnState_bkup);
	rande("state",randeState_bkup);
	randg("state",randgState_bkup);
	randp("state",randpState_bkup);
return;
end

%!test
%!	seedPrm = demoFunc0101_genSeedPrm("easy");
%!	seedPrm.sizeX = 100;
%!	seedPrm.sizeF = seedPrm.sizeX;
%!	funcPrm = demoFunc0101_genFuncPrm(seedPrm);
