	clear;
	commondefs;
	thisFile = "test_numoptCalcFullStep11";
	msg( thisFile, __LINE__, "Starting test." );
	setprngstates(0);
	tic();
	%
	msg( thisFile, __LINE__, "This is a hard 'has negative' divergent case." );
	%
	omega0 = 100;
	vecG = [ 2; 1 ];
	matH = [ 1, 2; 2, 1 ];
	msg( thisFile, __LINE__, sprintf( "omega0 = %e.", omega0 ) );
	%
	prm_calcFullStep = [];
	%%%prm_calcFullStep.deltaNormThresh = 1e-1;
	vecDelta = numoptCalcFullStep( omega0, vecG, matH, prm_calcFullStep );
	omega1 = omega0 + vecG'*vecDelta + 0.5*vecDelta'*matH*vecDelta;
	msg( thisFile, __LINE__, sprintf( "delta1Norm = %e;  omega1 = %e.", norm(vecDelta), omega1 ) );
	%assert( omega1 <= eps*omega0 )
	%
	coeffG = (vecG'*vecG)/(vecG'*matH*vecG);
	vecDeltaG = -coeffG*vecG;
	omegaG = omega0 + vecG'*vecDeltaG + 0.5*vecDeltaG'*matH*vecDeltaG;
	msg( thisFile, __LINE__, sprintf( "deltaGNorm = %e;  omegaG = %e.", norm(vecDeltaG), omegaG ) );
	%
	matI = eye(size(vecG,1));
	hFrobNorm = sqrt(sum(sum(matH.^2)));
	vecDeltaF = -( matH + hFrobNorm*matI ) \ vecG;
	omegaF = omega0 + vecG'*vecDeltaF + 0.5*vecDeltaF'*matH*vecDeltaF;
	msg( thisFile, __LINE__, sprintf( "deltaFNorm = %e;  omegaF = %e.", norm(vecDeltaF), omegaF ) );
	%
	msg( thisFile, __LINE__, "Finished test." );
	toc();
	%
	return
	
	numVals = 1000;
	muVals = 10.0.^linspace(0.0,log10(3.0*hFrobNorm),numVals);
	iotaVals = 0*muVals;
	iotaModelVals = 0*muVals;
	for n=1:numVals
		mu = muVals(n);
		[ matR, cholFlag ] = chol( matH + mu*matI );
		if (~cholFlag)
			vecDelta = -( matR \ (matR'\vecG) );
			omega = omega0 + vecDelta'*vecG + 0.5*vecDelta'*matH*vecDelta;
			iotaVals(n) = omega0-omega;
		end
		vecDeltaModel = -vecG/mu;
		omegaModel = omega0 + vecDeltaModel'*vecG + 0.5*vecDeltaModel'*matH*vecDeltaModel;
		iotaModelVals(n) = omega0-omegaModel;
	end
	%
	figure(1);
	loglog( muVals, iotaVals, 'o-' );
	grid on;
	%
	figure(2);
	plot( asinh(muVals), asinh(iotaVals-iotaModelVals), 'o-' );
	grid on;
	
	
