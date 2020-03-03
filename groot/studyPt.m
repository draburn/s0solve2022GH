function [ vecDelta_suggested, retCode, datOut ] = studyPt( ...
  funchF, vecX0, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	commondefs;
	thisFile = "studyPt";
	startTime = time();
	retCode = RETCODE__NOT_SET;
	%
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	reportInterval = mygetfield( prm, "reportInterval", 0.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	reportTimePrev = startTime - 0.1;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% SPECIFIC INIT.
	%
	sizeX = size(vecX0,1);
	assert( isrealarray(vecX0,[sizeX,1]) );
	vecF0 = funchF(vecX0);
	sizeF = size(vecF0,1);
	assert( isrealarray(vecF0,[sizeF,1]));
	funchJ = prm.funchJ;
	matJ0 = funchJ(vecX0);
	assert( isrealarray(matJ0,[sizeF,sizeX]) );
	%
	considerScl = mygetfield( prm, "considerScl", true );
	considerGradCurve = mygetfield( prm, "considerGradCurve", true );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PREP WORK.
	%
	vecG0 = matJ0'*vecF0;
	matH0 = matJ0'*matJ0;
	vecDiagH0 = diag(matH0);
	assert( 0.0 < min(vecDiagH0) );
	matD0 = diag(vecDiagH0);
	matI0 = eye(sizeX,sizeX);
	%
	%
	vecDelta_newton = -( matH0 \ vecG0 );
	funchDelta_newton = @(s)( vecDelta_newton * s );
	minScanPrm.funchDeltaSupportsMultiArg = true;
	minScanPrm.numPts1 = 21;
	s_newton_omegaMin = minScan( ...
	  funchF, vecX0, funchDelta_newton, minScanPrm );
	vecDelta_newton_omegaMin = funchDelta_newton(s_newton_omegaMin);
	clear minScanPrm;
	%
	%
	vecTemp = -vecG0;
	fTemp = vecTemp' * matH0 * vecTemp;
	assert( 0.0 < fTemp );
	vecDelta_gradDir = vecTemp * (-(vecTemp'*vecG0)/fTemp);
	funchDelta_gradDir = @(s)( vecDelta_gradDir * s );
	minScanPrm.funchDeltaSupportsMultiArg = true;
	minScanPrm.numPts1 = 21;
	s_gradDir_omegaMin = minScan( ...
	  funchF, vecX0, funchDelta_gradDir, minScanPrm );
	vecDelta_gradDir_omegaMin = funchDelta_gradDir(s_gradDir_omegaMin);
	clear minScanPrm;
	clear fTemp;
	clear vecTemp;
	%
	%
	if (considerScl)
		vecTemp = -( matD0 \ vecG0 );
		fTemp = vecTemp' * matH0 * vecTemp;
		assert( 0.0 < fTemp );
		vecDelta_gradDirScl = vecTemp * (-(vecTemp'*vecG0)/fTemp);
		funchDelta_gradDirScl = @(s)( vecDelta_gradDirScl * s );
		minScanPrm.funchDeltaSupportsMultiArg = true;
		minScanPrm.numPts1 = 21;
		s_gradDirScl_omegaMin = minScan( ...
		  funchF, vecX0, funchDelta_gradDirScl, minScanPrm );
		vecDelta_gradDirScl_omegaMin = funchDelta_gradDirScl(s_gradDirScl_omegaMin);
		clear minScanPrm
		clear fTemp;
		clear vecTemp;
	else
		vecDelta_gradDirScl = [];
		funcDelta_gradDirScl = -1.0;
		s_gradDirScl_omegaMin = [];
		vecDelta_gradDirScl_omegaMin = [];
	end
	%
	%
	mu0 = max(vecDiagH0);
	matA0 = matH0 + (mu0*matI0);
	matL0 = mu0 * matI0;
	funchDelta_levenberg = @(s)( -s * (( matA0 - (s*matL0))\vecG0) );
	minScanPrm.funchDeltaSupportsMultiArg = false;
	% But, could probably modify to make funchDeltaSupportsMultiArg true,
	% assuming we perform eig().
	minScanPrm.numPts1 = 201;
	s_levenberg_omegaMin = minScan( ...
	  funchF, vecX0, funchDelta_levenberg, minScanPrm );
	vecDelta_levenberg_omegaMin = funchDelta_levenberg(s_levenberg_omegaMin);
	clear minScanPrm;
	clear matL0;
	clear matA0;
	clear mu0;
	%
	%
	if (considerScl)
		matA0 = matH0 + matD0;
		matL0 = matD0;
		funchDelta_levenbergScl = @(s)( -s * (( matA0 - (s*matL0))\vecG0) );
		minScanPrm.funchDeltaSupportsMultiArg = false;
		% But, could probably modify to make funchDeltaSupportsMultiArg true,
		% assuming we perform eig().
		minScanPrm.numPts1 = 201;
		s_levenbergScl_omegaMin = minScan( ...
		  funchF, vecX0, funchDelta_levenbergScl, minScanPrm );
		vecDelta_levenbergScl_omegaMin = funchDelta_levenbergScl(s_levenbergScl_omegaMin);
		clear minScanPrm;
		clear matL0;
		clear matA0;
		clear mu0;
	else
		funchDelta_levenbergScl = [];
		s_levenbergScl_omegaMin = -1.0;
		vecDelta_levenbergScl_omegaMin = [];
	end
	%
	%
	if (considerGradCurve)
		[ matPsi, matLambda ] = eig( matH0 );
		assert( sum(sum(abs(((matPsi')*matPsi)-matI0))) < 10.0*(sizeX^3)*(eps^0.75) );
		assert( sum(sum(abs((matPsi*(matPsi'))-matI0))) < 10.0*(sizeX^3)*(eps^0.75) );
		% Is (sizeX^3) correct?
		vecPsiTN = matPsi'*(-matH0\vecG0);
		lambdaMin = min(diag(matLambda));
		matSigma = matLambda / lambdaMin;
		funchDelta_gradCurve = @(s)( ...
		  matPsi * ( vecPsiTN - (diag((1.0-s).^diag(matSigma))*vecPsiTN) ) );
		minScanPrm.funchDeltaSupportsMultiArg = false;
		% But, could probably modify to make funchDeltaSupportsMultiArg true.
		minScanPrm.numPts1 = 201;
		s_gradCurve_omegaMin = minScan( ...
		  funchF, vecX0, funchDelta_gradCurve, minScanPrm );
		vecDelta_gradCurve_omegaMin = funchDelta_gradCurve(s_gradCurve_omegaMin);
		clear matSigma;
		clear lambdaMin;
		clear vecPsiTN;
		clear matLambda;
		clear matPsi;
	else
		funchDelta_gradCurve = [];
		s_gradCurve_omgeaMin = -1.0;
		vecDelta_gradCurve_omegaMin = [];
	end
	%
	%
	if (considerGradCurve&&considerScl)
		matD0InvSqrt = diag(1./sqrt(vecDiagH0));
		matH0Scl = matD0InvSqrt * matH0 * matD0InvSqrt;
		vecG0Scl = matD0InvSqrt * vecG0;
		[ matPsi, matLambda ] = eig( matH0Scl );
		assert( sum(sum(abs(((matPsi')*matPsi)-matI0))) < 10.0*(sizeX^3)*(eps^0.75) );
		assert( sum(sum(abs((matPsi*(matPsi'))-matI0))) < 10.0*(sizeX^3)*(eps^0.75) );
		% Is (sizeX^3) correct?
		vecPsiTN = -matPsi'*(matH0Scl\vecG0Scl);
		vecDiagLambda = diag(matLambda);
		lambdaMin = min(vecDiagLambda);
		assert( 0.0 < lambdaMin );
		vecDiagSigma = vecDiagLambda / lambdaMin;
		funchDelta_gradCurveScl = @(s)( ...
		  vecDelta_newton - (matPsi*( ((1.0-s).^vecDiagSigma) .* vecPsiTN))  );
		minScanPrm.funchDeltaSupportsMultiArg = false;
		% But, could probably modify to make funchDeltaSupportsMultiArg true.
		minScanPrm.numPts1 = 201;
		s_gradCurveScl_omegaMin = minScan( ...
		  funchF, vecX0, funchDelta_gradCurveScl, minScanPrm );
		vecDelta_gradCurveScl_omegaMin = funchDelta_gradCurveScl(s_gradCurveScl_omegaMin);
		clear vecDiagSigma;
		clear lambdaMin;
		clear vecDiagLambda;
		clear vecPsiTN;
		clear matLambda;
		clear matPsi;
		clear vecG0Scl;
		clear matH0Scl;
		clear matD0InvSqrt;
	else
		funcDelta_gradCirveScl = [];
		s_gradCurveScl_omegaMin = -1.0;
		vecDelta_gradCurveScl_omegaMin = [];
	end
	%
	%
	vecDelta_suggested = vecDelta_levenberg_omegaMin;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% POPULATE DATOUT.
	%
	datOut.funchF = funchF;
	datOut.vecX0 = vecX0;
	datOut.prm = prm;
	datOut.startTime = startTime;
	%
	datOut.verbLev = verbLev;
	datOut.reportInterval = reportInterval;
	%
	datOut.sizeX = sizeX;
	datOut.sizeF = sizeF;
	datOut.vecF0 = vecF0;
	datOut.matJ0 = matJ0;
	%
	datOut.considerScl = considerScl;
	datOut.considerGradCurve = considerGradCurve;
	%
	datOut.vecG0 = vecG0;
	datOut.matH0 = matH0;
	datOut.vecDiagH0 = vecDiagH0;
	datOut.matD0 = matD0;
	datOut.matI0 = matI0;
	%
	datOut.vecDelta_newton = vecDelta_newton;
	datOut.vecDelta_gradDir = vecDelta_gradDir;
	datOut.vecDelta_gradDirScl = vecDelta_gradDirScl;
	%
	datOut.funchDelta_newton = funchDelta_newton;
	datOut.funchDelta_gradDir = funchDelta_gradDir;
	datOut.funchDelta_gradDirScl = funchDelta_gradDirScl;
	datOut.funchDelta_levenberg = funchDelta_levenberg;
	datOut.funchDelta_levenbergScl = funchDelta_levenbergScl;
	datOut.funchDelta_gradCurve = funchDelta_gradCurve;
	datOut.funchDelta_gradCurveScl = funchDelta_gradCurveScl;
	%
	datOut.vecDelta_newton_omegaMin = vecDelta_newton_omegaMin;
	datOut.vecDelta_gradDir_omegaMin = vecDelta_gradDir_omegaMin;
	datOut.vecDelta_gradDirScl_omegaMin = vecDelta_gradDirScl_omegaMin;
	datOut.vecDelta_levenberg_omegaMin = vecDelta_levenberg_omegaMin;
	datOut.vecDelta_levenbergScl_omegaMin = vecDelta_levenbergScl_omegaMin;
	datOut.vecDelta_gradCurve_omegaMin = vecDelta_gradCurve_omegaMin;
	datOut.vecDelta_gradCurveScl_omegaMin = vecDelta_gradCurveScl_omegaMin;
	%
	datOut.s_newton_omegaMin = s_newton_omegaMin;
	datOut.s_gradDir_omegaMin = s_gradDir_omegaMin;
	datOut.s_gradDirScl_omegaMin = s_gradDirScl_omegaMin;
	datOut.s_levenberg_omegaMin = s_levenberg_omegaMin;
	datOut.s_levenbergScl_omegaMin = s_levenbergScl_omegaMin;
	datOut.s_gradCurve_omegaMin = s_gradCurve_omegaMin;
	datOut.s_gradCurveScl_omegaMin = s_gradCurveScl_omegaMin;
	%
	msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
	 "Norm delta values: %g, %g, %g; %g, %g, %g, %g, %g, %g, %g.", ...
	 norm(vecDelta_newton), ...
	 norm(vecDelta_gradDir), ...
	 norm(vecDelta_gradDirScl), ...
	 norm(vecDelta_newton_omegaMin), ...
	 norm(vecDelta_gradDir_omegaMin), ...
	 norm(vecDelta_gradDirScl_omegaMin), ...
	 norm(vecDelta_levenberg_omegaMin), ...
	 norm(vecDelta_levenbergScl_omegaMin), ...
	 norm(vecDelta_gradCurve_omegaMin), ...
	 norm(vecDelta_gradCurveScl_omegaMin) ));
	msg_progress( verbLev, thisFile, __LINE__, sprintf( ...
	 "S values: %g, %g, %g, %g, %g, %g, %g.", ...
	 s_newton_omegaMin, ...
	 s_gradDir_omegaMin, ...
	 s_gradDirScl_omegaMin, ...
	 s_levenberg_omegaMin, ...
	 s_levenbergScl_omegaMin, ...
	 s_gradCurve_omegaMin, ...
	 s_gradCurveScl_omegaMin ));
	%
return;
end

%!test
%!	test_studyPt;
