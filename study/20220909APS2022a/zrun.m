clear
printf("\n\n");
mydefs;
startTime = time();
%
seedTime = mod( round(now*1E11), 1E8 ); % In case you want this.
%
probSetPrm = [];
probSetPrm.numProbs = 10;
probSetPrm.probType = "test1";
probSetPrm.numUnknowns = 20;
probSetPrm.setSeed = 0;
%
%
solverPrm = [];
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
%solverPrm.solverFunch = @groot_fsolve;
solverPrm.solverFunch = @groot_jfnk_basic;
%solverPrm.solverFunch = @groot_jfnk_basic; solverPrm.btCoeff = 0.0;
%
%solverPrm.btCoeff = 0.1;
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
	[ vecXBest, strGrootFlag, fevalCount, solverDatOut ] = solverPrm.solverFunch( funchF, vecX0, solverPrm );
	%
	assert( isrealarray(vecXBest,[probDat.sizeX,1]) );
	assert( ischar(strGrootFlag) );
	assert( isvector(strGrootFlag) );
	assert( isposintscalar(fevalCount) );
	%
	vecFBest = funchF(vecXBest);
	assert( isrealarray(vecFBest,[probDat.sizeF,1]) );
	fBest = norm(vecFBest);
	switch ( strGrootFlag )
	case STR_GROOT_FLAG__CNVG
		assert( fBest <= solverPrm.fTol*(1.0+100.0*eps) );
		assert( fevalCount <= solverPrm.fevalLimit );
		zrunDat.succCount++;
		zrunDat.succFevalVals = [ zrunDat.succFevalVals, fevalCount ];
	case { STR_GROOT_FLAG__STOP, STR_GROOT_FLAG__STALL }
		assert( fBest >= solverPrm.fTol*(1.0-100.0*eps) );
		zrunDat.failCount++;
	otherwise
		error(["Unsupported value of strGrootFlag (\"" strGrootFlag "\")."]);
	endswitch
	msg( __FILE__, __LINE__, sprintf("    %3d  %11s:   %7s  %6d  %10.3e", probIndex, probSizeStr, strGrootFlag, fevalCount, fBest ) );
	%
	zrunDat.prob(probIndex).vecXBest = vecXBest;
	zrunDat.prob(probIndex).fevalCount = fevalCount;
	zrunDat.prob(probIndex).matRecordX = solverDatOut.matRecordX;
	zrunDat.prob(probIndex).matInfoA = solverDatOut.matInfoA;
	zrunDat.prob(probIndex).matInfoB = mygetfield( solverDatOut, "matInfoB", [] );
	if ( max(size(probList)) == 1 )
		semilogy( solverDatOut.matInfoA(:,2), solverDatOutmatInfoA(:,4), 'o-') ;
		grid on;
	endif
endfor
%
zrun__finish;
%
msg( __FILE__, __LINE__, "Goodbye." );
printf("\n\n");
