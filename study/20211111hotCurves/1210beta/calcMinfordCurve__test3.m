clear;
thisFile = "calcMinfordCurve__test3";
commondefs;
setprngstates(0);
numFigs = 0;
startTime = time();
%

numTrials = 10;
msg( thisFile, __LINE__, sprintf( "Performing %d trials...", numTrials ) );
time0 = time();
for t=1:numTrials
	sizeX = 2+abs(round(10.0*randn())); % gSurf will be zero if sizeX = 1.
	omega0 = abs(randn());
	vecG0 = randn(sizeX,1);
	matHGen = randn(sizeX,sizeX);
	matH0 = matHGen'+matHGen;
	funchOmega = @(x)( omega0 + x'*vecG0 + 0.5*(x'*matH0*x) );
	funchG = @(x)( vecG0 + matH0*x );
	funchH = @(x)( matH0 );
	clear matHGen;
	%
	testPrm = [];
	if ( 1==t )
		msg( thisFile, __LINE__, sprintf( "Showing values for trial %d...", t ) );
		testPrm.verbLev = VERBLEV__PROGRESS;
		echo__omega0 = omega0
		echo__vecG0 = vecG0
		echo__matH0 = matH0
	end
	calcMinfordCurve__testCalc( sizeX, funchOmega, funchG, funchH, testPrm );
	msg( thisFile, __LINE__, sprintf( "Passed trial %d...", t ) );
end
msg( thisFile, __LINE__, sprintf( "Passed %d trials in %0.3fs.", numTrials, time()-time0 ) );
