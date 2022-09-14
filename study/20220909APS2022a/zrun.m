clear
printf("\n\n");
mydefs;
startTime = time();
%
seedTime = mod( round(now*1E11), 1E8 ); % In case you want this.
%
probSetPrm = [];
probSetPrm.numProbs = 20;
probSetPrm.probType = "test0";
probSetPrm.numUnknowns = 20;
probSetPrm.setSeed = 0;
%
solverPrm = [];
solverPrm.verbLev = VERBLEV__FLAGGED;
solverPrm.valdLev = VALDLEV__HIGH;
solverPrm.solverName = "fsolve";
solverPrm.iterLimit = 100;
solverPrm.fevalLimit = 100*probSetPrm.numUnknowns;
solverPrm.fTol = 1.0e-8;
solverPrm.fallTol = solverPrm.fTol / 100.0;
solverPrm.stepTol = 1.0e-8;
solverPrm.epsFD = eps^(1.0/3.0);
%
zrun__start;
%
for probIndex = 1 : probSetPrm.numProbs
	probGenPrm = [];
	[ funchF, vecX0, probDat ] = zrun_genProb( probSetPrm.probType, probSetPrm.numUnknowns, probSetPrm.probSeeds(probIndex), probGenPrm );
	assert( isrealarray(vecX0,[probDat.sizeX,1]) );
	probSizeStr = sprintf( "(%dx%d)", probDat.sizeF, probDat.sizeX );
	%
	[ vecXBest, fevalCount ] = groot_fsolve( funchF, vecX0, solverPrm );
	%[ vecXBest, fevalCount ] = groot_jfnk_crude( funchF, vecX0, solverPrm );
	assert( isrealarray(vecXBest,[probDat.sizeX,1]) );
	assert( isposintscalar(fevalCount) );
	vecFBest = funchF(vecXBest);
	assert( isrealarray(vecFBest,[probDat.sizeF,1]) );
	fBest = norm(vecFBest);
	if ( fBest <= solverPrm.fTol )
		resultStr = " success ";
	else
		resultStr = "*FAILURE*";
	endif
	msg( __FILE__, __LINE__, sprintf("    %3d  %11s:   %9s  %6d  %10.3e", probIndex, probSizeStr, resultStr, fevalCount, fBest ) );
endfor
%
zrun__finish;
%
msg( __FILE__, __LINE__, "Goodbye." );
printf("\n\n");
