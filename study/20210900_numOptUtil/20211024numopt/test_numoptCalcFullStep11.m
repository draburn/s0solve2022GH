	clear;
	commondefs;
	thisFile = "test_numoptCalcFullStep11";
	msg( thisFile, __LINE__, "Starting test." );
	setprngstates(0);
	tic();
	%
	msg( thisFile, __LINE__, "This is an easy 'has negative' divergent case." );
	%
	omega0 = 1;
	vecG = [ 2; 1 ];
	matH = [ 1, 2; 2, 1 ];
	msg( thisFile, __LINE__, sprintf( "omega0 = %e.", omega0 ) );
	%
	vecDelta = numoptCalcFullStep( omega0, vecG, matH );
	omega1 = omega0 + vecG'*vecDelta + 0.5*vecDelta'*matH*vecDelta;
	msg( thisFile, __LINE__, sprintf( "deltaNorm = %e;  omega1 = %e.", norm(vecDelta), omega1 ) );
	assert( omega1 <= eps*omega0 )
	%
	msg( thisFile, __LINE__, "Finished test." );
	toc();
