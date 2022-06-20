if (isempty(r.prmMemo))
	r.runName = sprintf( "%d", r.runType );
else
	r.runName = sprintf( "%d %s", r.runType, r.prmMemo );
endif
msg( __FILE__, __LINE__, sprintf( "Starting run %d/%d: '%s'...", runIndex, numRuns, r.runName ) );
runStartTime = time();
%
switch( r.runType )
case 50
	r.runTypeDescrip = "Quess N-R + minscan";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_050( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 444
	r.runTypeDescrip = "444";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_444( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 550
	r.runTypeDescrip = "Quess JFNK + minscan";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_550( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 800
	r.runTypeDescrip = "JFNK+TR+AP";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_800( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 802
	r.runTypeDescrip = "800tweaky";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_800tweaky( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 805
	r.runTypeDescrip = "800simple";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_800simple( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 810
	r.runTypeDescrip = "800splitspace";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_800splitspace( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 820
	r.runTypeDescrip = "800scomb";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_800scomb( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 85551
	r.runTypeDescrip = "800sssl";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_800sssl( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 815353
	r.runTypeDescrip = "800insensed";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_800insensed( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 8153532
	r.runTypeDescrip = "800insensed2";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_800insensed2( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 904
	r.runTypeDescrip = "phiPatch?";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_904( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 940
	r.runTypeDescrip = "phiPatch + slinsolf?";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_940( vecX0, funchF, r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.iterCount;
	r.fevalCountOfStep = r.datOut.fevalCountVals;
	r.fBestNormOfStep = r.datOut.fNormVals;
	r.isValid = true;
case 1100
	r.runTypeDescrip = "z100";
	[ r.vecXF, r.vecFF, r.datOut ] = zlinsolf100( funchF, vecX0, [], r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.stepCount;
	r.fevalCountOfStep = r.datOut.fevalCountOfSteps;
	r.fBestNormOfStep = r.datOut.fNormOfSteps;
	%%%r.fevalCountOfStep = r.datOut.fevalCountVals(r.datOut.iterCountOfSteps);
	%%%r.fBestNormOfStep = r.datOut.fNormVals(r.datOut.iterCountOfSteps);
	r.isValid = true;
case 1150
	r.runTypeDescrip = "z150";
	[ r.vecXF, r.vecFF, r.datOut ] = zlinsolf150( funchF, vecX0, [], r.prm );
	r.fevalCount = r.datOut.fevalCount;
	r.stepCount = r.datOut.stepCount;
	r.fevalCountOfStep = r.datOut.fevalCountOfSteps;
	r.fBestNormOfStep = r.datOut.fNormOfSteps;
	r.isValid = true;
case 1195
	r.runTypeDescrip = "z195";
	[ r.vecXF, r.vecFF, r.retCode, r.fevalCount, r.stepsCount, r.datOut ] = zlinsolf195( funchF, vecX0, [], r.prm );
	r.stepCount = r.stepsCount;
	r.fevalCountOfStep = r.datOut.fevalCountOfSteps;
	r.fBestNormOfStep = r.datOut.fNormOfSteps;
	r.isValid = true;
case 2100
	r.runTypeDescrip = "s100";
	[ r.vecXF, r.vecFF, r.retCode, r.fevalCount, r.stepsCount, r.datOut ] = sxsolf100( funchF, vecX0, [], r.prm );
	r.stepCount = r.stepsCount;
	r.fevalCountOfStep = r.datOut.fevalCountOfSteps;
	r.fBestNormOfStep = r.datOut.fNormOfSteps;
	r.isValid = true;
case 2114
	r.runTypeDescrip = "s114";
	[ r.vecXF, r.vecFF, r.retCode, r.fevalCount, r.stepsCount, r.datOut ] = sxsolf114( funchF, vecX0, [], r.prm );
	r.stepCount = r.stepsCount;
	r.fevalCountOfStep = r.datOut.fevalCountOfSteps;
	r.fBestNormOfStep = r.datOut.fNormOfSteps;
	r.isValid = true;
case 2124
	r.runTypeDescrip = "s124";
	[ r.vecXF, r.vecFF, r.retCode, r.fevalCount, r.stepsCount, r.datOut ] = sxsolf124( funchF, vecX0, [], r.prm );
	r.stepCount = r.stepsCount;
	r.fevalCountOfStep = r.datOut.fevalCountOfSteps;
	r.fBestNormOfStep = r.datOut.fNormOfSteps;
	r.isValid = true;
case 2180
	r.runTypeDescrip = "s180";
	[ r.vecXF, r.vecFF, r.retCode, r.fevalCount, r.stepsCount, r.datOut ] = sxsolf180( funchF, vecX0, [], r.prm );
	r.stepCount = r.stepsCount;
	r.fevalCountOfStep = r.datOut.fevalCountOfSteps;
	r.fBestNormOfStep = r.datOut.fNormOfSteps;
	r.isValid = true;
case 2181
	r.runTypeDescrip = "s181";
	[ r.vecXF, r.vecFF, r.retCode, r.fevalCount, r.stepsCount, r.datOut ] = sxsolf181( funchF, vecX0, [], r.prm );
	r.stepCount = r.stepsCount;
	r.fevalCountOfStep = r.datOut.fevalCountOfSteps;
	r.fBestNormOfStep = r.datOut.fNormOfSteps;
	r.isValid = true;
case 2182
	r.runTypeDescrip = "s182";
	[ r.vecXF, r.vecFF, r.retCode, r.fevalCount, r.stepsCount, r.datOut ] = sxsolf182( funchF, vecX0, [], r.prm );
	r.stepCount = r.stepsCount;
	r.fevalCountOfStep = r.datOut.fevalCountOfSteps;
	r.fBestNormOfStep = r.datOut.fNormOfSteps;
	r.isValid = true;
case 2183
	r.runTypeDescrip = "s183";
	[ r.vecXF, r.vecFF, r.retCode, r.fevalCount, r.stepsCount, r.datOut ] = sxsolf183( funchF, vecX0, [], r.prm );
	r.stepCount = r.stepsCount;
	r.fevalCountOfStep = r.datOut.fevalCountOfSteps;
	r.fBestNormOfStep = r.datOut.fNormOfSteps;
	r.isValid = true;
case 21801
	r.runTypeDescrip = "s181wsparsep";
	[ r.vecXF, r.vecFF, r.retCode, r.fevalCount, r.stepsCount, r.datOut ] = sxsolf181wsparsep( funchF, vecX0, [], r.prm );
	r.stepCount = r.stepsCount;
	r.fevalCountOfStep = r.datOut.fevalCountOfSteps;
	r.fBestNormOfStep = r.datOut.fNormOfSteps;
	r.isValid = true;
case 21802
	r.runTypeDescrip = "s181perCompareAP";
	[ r.vecXF, r.vecFF, r.retCode, r.fevalCount, r.stepsCount, r.datOut ] = sxsolf181perCompareAP( funchF, vecX0, [], r.prm );
	r.stepCount = r.stepsCount;
	r.fevalCountOfStep = r.datOut.fevalCountOfSteps;
	r.fBestNormOfStep = r.datOut.fNormOfSteps;
	r.isValid = true;
otherwise
	msg( __FILE__, __LINE__, sprintf( "ERROR: Invalid runType (%d)", r.runType ) );
	r.isValid = false;
	return;
endswitch
markerTypes = "+o*xsd^v<>ph";
r.mlStyle = [ markerTypes(mod(runIndex,length(markerTypes))+1), "-" ];
r.mSize = 10+3*(numRuns-runIndex);
r.stepCountOfStep = ( 0 : r.stepCount );
%
r.elapsedTime = time()-runStartTime;
msg( __FILE__, __LINE__, sprintf( "Run '%s' (%s) reached %0.3e in %0.3es with %d fevals.", ...
  r.runName, r.runTypeDescrip, norm(r.vecFF), r.elapsedTime, r.fevalCount ) );
