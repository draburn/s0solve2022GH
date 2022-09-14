clear
printf("\n\n");
mydefs;
startTime = time();
%
seedTime = mod( round(now*1E11), 1E8 ); % In case you want this.
%
probSetPrm = [];
probSetPrm.numProbs = 10;
probSetPrm.probType = "testr";
probSetPrm.numUnknowns = 10;
probSetPrm.setSeed = 0;
%
%
solverPrm = [];
%solverPrm.solverFunch = @groot_fsolve;
solverPrm.solverFunch = @groot_jfnk_crude;
%
solverPrm.verbLev = VERBLEV__FLAGGED;
solverPrm.valdLev = VALDLEV__HIGH;
solverPrm.iterLimit = 100;
solverPrm.fevalLimit = 100*probSetPrm.numUnknowns;
solverPrm.fTol = 1.0e-8;
solverPrm.fallTol = solverPrm.fTol / 100.0;
solverPrm.stepTol = 1.0e-8;
solverPrm.epsFD = eps^(1.0/3.0);

%
%solverPrm.linsolfPrm.tol = 1.0e-2;
%
%
zrun__start;
%
probList = ( 1 : probSetPrm.numProbs );
%probList = [3];
for probIndex = probList
	if ( stopsignalpresent() )
		msg( __FILE__, __LINE__, "Received stopsignal." );
		return;
	endif
	%
	probGenPrm = [];
	[ funchF, vecX0, probDat ] = zrun_genProb( probSetPrm.probType, probSetPrm.numUnknowns, probSetPrm.probSeeds(probIndex), probGenPrm );
	assert( isrealarray(vecX0,[probDat.sizeX,1]) );
	probSizeStr = sprintf( "(%dx%d)", probDat.sizeF, probDat.sizeX );
	%
	[ vecXBest, fevalCount, matCnvg ] = solverPrm.solverFunch( funchF, vecX0, solverPrm );
	%
	assert( isrealarray(vecXBest,[probDat.sizeX,1]) );
	assert( isposintscalar(fevalCount) );
	vecFBest = funchF(vecXBest);
	assert( isrealarray(vecFBest,[probDat.sizeF,1]) );
	fBest = norm(vecFBest);
	if ( fBest <= solverPrm.fTol )
		resultStr = " success ";
		zrunDat.succCount++;
		zrunDat.succFevalVals = [ zrunDat.succFevalVals, fevalCount ];
	else
		resultStr = "*FAILURE*";
		zrunDat.failCount++;
	endif
	msg( __FILE__, __LINE__, sprintf("    %3d  %11s:   %9s  %6d  %10.3e", probIndex, probSizeStr, resultStr, fevalCount, fBest ) );
	%
	zrunDat.prob(probIndex).vecXBest = vecXBest;
	zrunDat.prob(probIndex).fevalCount = fevalCount;
	zrunDat.prob(probIndex).matCnvg = matCnvg;
	%
	if ( max(size(probList)) == 1 )
		semilogy( matCnvg(:,1), matCnvg(:,2), 'o-') ;
		grid on;
	endif
endfor
%
zrun__finish;
%
msg( __FILE__, __LINE__, "Goodbye." );
printf("\n\n");
