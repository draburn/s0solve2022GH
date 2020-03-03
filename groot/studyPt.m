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
	matD0 = diag(vecDiagH0);
	matI0 = eye(sizeX,sizeX);
	%
	%
	vecDelta_newton = -( matH0 \ vecG0 );
	funchDelta_newton = @(s)( vecDelta_newton * s );
	minScanPrm.numPts1 = 21;
	minScanPrm.funchDeltaSupportsMultiArg = true;
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
	minScanPrm.numPts1 = 21;
	minScanPrm.funchDeltaSupportsMultiArg = true;
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
		minScanPrm.numPts1 = 21;
		minScanPrm.funchDeltaSupportsMultiArg = true;
		s_gradDirScl_omegaMin = minScan( ...
		  funchF, vecX0, funchDelta_gradDirScl, minScanPrm );
		vecDelta_gradDirScl_omegaMin = funchDelta_gradDirScl(s_gradDirScl_omegaMin);
		clear minScanPrm
		clear fTemp;
		clear vecTemp;
	else
		vecDelta_gradDirScl = [];
		funcDelta_gradDirScl = [];
		s_gradDirScl_omegaMin = [];
		vecDelta_gradDirScl_omegaMin = [];
	end
	%
	%
	%
	%funchDeltaLevenberg = @(s)( -s*( matH0 + 
	%vecXLevenbergLin = vecXNewtonLin;  % FWIW;
	%
	%
	%
	if (considerScl)
		%funchDeltaLevenbergScl = @(s)(...);
		%vecXLevenbergSclLin = vecXNewtonLin; % FWIW;
	end
	%
	%
	%
	if (considerGradCurve)
		%funchDeltaGradCurve = @(s)( ... );
		%vecXGradCurveLin = vecXNewtonLin;  % FWIW;
	end
	%
	%
	%
	if (considerGradCurve&&considerScl)
		%funchDeltaGradCurveScl = @(s)( ... );
		%vecXGradCurveSclLin = vecXNewtonLin; % FWIW;
	end
	%
	%
	%
	%
	%vecXGradDirOmega = ...
	%
	if (considerScl)
		%vecXGradDirSclOmega = ...
	end
	%
	%vecXLevenbergOmega = ...
	%
	if (considerScl)
		%vecXLevenbergSclOmega = ...
	end
	%
	if (considerGradCurve)
		%vecXGradCurveOmega = ...
	end
	%
	if (considerGradCurve&&considerScl)
		%vecXGradCurveSclOmega = ...
	end
	%
	%
	vecDelta_suggested = vecDelta_newton_omegaMin;
	%
	%
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
	datOut.vecDelta_newton_omegaMin = vecDelta_newton_omegaMin;
	%datOut.vecDelta_levenberg_omegaMin = vecDelta_levenberg_omegaMin;
	%datOut.vecDelta_levenbergScl_omegaMin = vecDelta_levenbergScl_omegaMin;
	datOut.vecDelta_gradDir_omegaMin = vecDelta_gradDir_omegaMin;
	datOut.vecDelta_gradDirScl_omegaMin = vecDelta_gradDirScl_omegaMin;
	%datOut.vecDelta_gradCurve_omegaMin = vecDelta_gradCurve_omegaMin;
	%datOut.vecDelta_gradCurveScl_omegaMin = vecDelta_gradCurveScl_omegaMin;
	%
	datOut.s_newton_omegaMin = s_newton_omegaMin;
	%datOut.s_levenberg_omegaMin = s_levenberg_omegaMin;
	%datOut.s_levenbergScl_omegaMin = s_levenbergScl_omegaMin;
	datOut.s_gradDir_omegaMin = s_gradDir_omegaMin;
	datOut.s_gradDirScl_omegaMin = s_gradDirScl_omegaMin;
	%datOut.s_gradCurve_omegaMin = s_gradCurve_omegaMin;
	%datOut.s_gradCurveScl_omegaMin = s_gradCurveScl_omegaMin;
return;
end

%!test
%!	test_studyPt;
