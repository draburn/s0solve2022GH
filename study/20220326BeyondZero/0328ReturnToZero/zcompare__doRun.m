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
	msg( __FILE__, __LINE__, sprintf( "ERROR: Invalid runType (%d)", runType ) );
	runDat = runDat_dummy;
	return;
endswitch
%
r.elapsedTime = time()-runStartTime;
msg( __FILE__, __LINE__, sprintf( "'%s' reached %0.3e in %0.3es with %d fevals.", ...
  r.runTypeDescrip, norm(r.vecFF), r.elapsedTime, r.datOut.fevalCount ) );
