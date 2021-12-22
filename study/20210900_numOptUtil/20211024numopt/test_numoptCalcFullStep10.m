	clear;
	commondefs;
	thisFile = "test_numoptCalcFullStep4";
	msg( thisFile, __LINE__, "Starting test." );
	setprngstates(0);
	tic();
	%
	msg( thisFile, __LINE__, "This is a positive semi-definite divergent case." );
	%
	omega0 = 10;
	vecG = [ 1; 1 ];
	matH = [ 1, 0; 0, 0 ];
	msg( thisFile, __LINE__, sprintf( "omega0 = %e.", omega0 ) );
	%
	vecDelta = numoptCalcFullStep( omega0, vecG, matH );
	omega1 = omega0 + vecG'*vecDelta + 0.5*vecDelta'*matH*vecDelta;
	msg( thisFile, __LINE__, sprintf( "omega1 = %e.", omega1 ) );
	assert( omega1 <= eps*omega0 )
	%
	msg( thisFile, __LINE__, "Finished test." );
	toc();
