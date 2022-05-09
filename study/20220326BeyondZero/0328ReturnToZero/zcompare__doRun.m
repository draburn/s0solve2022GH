runStartTime = time();
%
switch( r.runType )
case 50
	r.runTypeDescrip = "Quess N-R + minscan";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_050( vecX0, funchF, r.prm );
case 550
	r.runTypeDescrip = "Quess JFNK + minscan";
	[ r.vecXF, r.vecFF, r.datOut ] = findZero_550( vecX0, funchF, r.prm );
case 1100
	r.runTypeDescrip = "z100";
	[ r.vecXF, r.vecFF, r.datOut ] = zlinsolf100( funchF, vecX0, [], r.prm );
otherwise
	msg( __FILE__, __LINE__, sprintf( "ERROR: Invalid runType (%d)", r.runType ) );
	r.vecFF = funchF(vecX0);
	r.runTypeDescrip = "INVALID";
	r.elapsedTime = 0.0;
	r.stepCount = -1;
	r.fevalCount = -1;
	return;
endswitch
%
r.elapsedTime = time()-runStartTime;
if (isempty(mygetfield(r.datOut,"stepCount",[])))
	r.stepCount = r.datOut.iterCount;
else
	r.stepCount = r.datOut.stepCount;
endif
r.fevalCount = r.datOut.fevalCount;
msg( __FILE__, __LINE__, sprintf( "Run '%s' (%s) reached %0.3e in %0.3es with %d fevals.", ...
  r.runName, r.runTypeDescrip, norm(r.vecFF), r.elapsedTime, r.datOut.fevalCount ) );
