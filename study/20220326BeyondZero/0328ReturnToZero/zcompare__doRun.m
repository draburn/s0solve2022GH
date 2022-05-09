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
case 550
	r.runTypeDescrip = "Quess JFNK + minscan";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_550( vecX0, funchF, r.prm );
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
	r.fevalCountOfStep = r.datOut.fevalCountVals(r.datOut.iterCountOfSteps);
	r.fBestNormOfStep = r.datOut.fNormVals(r.datOut.iterCountOfSteps);
	r.isValid = true;
otherwise
	msg( __FILE__, __LINE__, sprintf( "ERROR: Invalid runType (%d)", r.runType ) );
	r.isValid = false;
	return;
endswitch
%
r.elapsedTime = time()-runStartTime;
msg( __FILE__, __LINE__, sprintf( "Run '%s' (%s) reached %0.3e in %0.3es with %d fevals.", ...
  r.runName, r.runTypeDescrip, norm(r.vecFF), r.elapsedTime, r.datOut.fevalCount ) );
