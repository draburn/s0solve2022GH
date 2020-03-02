function [ vecXSuggested, retCode, datOut ] = studyPt( ...
  funchF, vecX0, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% COMMON INIT.
	%
	datOut.funchF = funchF;
	datOut.vecX0 = vecX0;
	datOut.prm = prm;
	%
	commondefs;
	thisFile = "studyPt";
	startTime = time();
	%
	datOut.startTime = startTime;
	%
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	reportInterval = mygetfield( prm, "reportInterval", 0.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	%
	datOut.verbLev = verbLev;
	datOut.reportInterval = reportInterval;
	%
	reportTimePrev = startTime - 0.1;
	%
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
	datOut.sizeX = sizeX;
	datOut.sizeF = sizeF;
	datOut.vecX0 = vecX0;
	datOut.vecF0 = vecF0;
	datOut.matJ0 = matJ0;
	%
	%
	%
	vecXSecret = prm.vecXSecret;
	assert( isrealarray(vecXSecret,[sizeX,1]) );
	%
	datOut.vecXSecret = vecXSecret;
	%
	%
	%
	considerPicard = mygetfield( prm, "considerPicard", true );
	considerScl = mygetfield( prm, "considerScl", true );
	considerGradCurve = mygetfield( prm, "considerGradCurve", true );
	%
	datOut.considerPicard = considerPicard;
	datOut.considerScl = considerScl;
	datOut.considerGradCurve = considerGradCurve;
	%
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% PREP WORK.
	%
	vecG0 = matJ0'*vecF0;
	matH0 = matJ0'*matJ0;
	vecD0 = diag(matH0);
	matD0 = diag(vecD0);
	matI0 = eye(sizeX,sizeX);
	matXF = eye(sizeF,sizeX);
	%
	datOut.vecG0 = vecG0;
	datOut.matH0 = matH0;
	datOut.vecD0 = vecD0; % Drop this?
	datOut.matD0 = matD0; % Drop this?
	datOut.matI0 = matI0; % Drop this?
	datOut.matXF = matXF; % Drop this?
	%
	%
	%
	vecDeltaNewton = -( matH0 \ vecG0 );
	vecXNewtonLin  = vecX0 + vecDeltaNewton;
	%vecDeltaNewtonScl = vecDeltaNewton; % FWIW
	%vecXNewtonSclLin  = vecXNewtonLin; % FWIW;
	%
	%vecXNewtonOmega = ...
	%
	datOut.vecDeltaNewton = vecDeltaNewton;
	datOut.vecXNewtonLin = vecXNewtonLin;
	%
	%
	%
	vecTemp = -vecG0;
	fTemp = vecTemp' * matH0 * vecTemp;
	assert( 0.0 < fTemp );
	vecDeltaGradDir = vecTemp * (-(vecTemp'*vecG0)/fTemp);
	vecXGradDirLin = vecX0 + vecDeltaGradDir;
	%
	%
	%
	datOut.vecDeltaGradDir = vecDeltaGradDir;
	datOut.vecXGradDirLin = vecXGradDirLin;
	%
	%
	%
	if (considerScl)
		vecTemp = -( matD0 \ vecG0 );
		fTemp = vecTemp' * matH0 * vecTemp;
		assert( 0.0 < fTemp );
		vecDeltaGradDirScl = vecTemp * (-(vecTemp'*vecG0)/fTemp);
		vecXGradDirSclLin = vecX0 + vecDeltaGradDirScl;
		%
		datOut.vecDeltaGradDirScl = vecDeltaGradDirScl;
		datOut.vecXGradDirSclLin = vecXGradDirSclLin;
	end
	%
	if (considerPicard)
		vecTemp = -( matXF * vecF0 );
		fTemp = vecTemp' * matH0 * vecTemp;
		assert( 0.0 < fTemp );
		vecDeltaPicard = vecTemp * (-(vecTemp'*vecG0)/fTemp);
		vecXPicardLin  = vecX0 + vecDeltaPicard;
		%
		datOut.vecDeltaPicard = vecDeltaPicard;
		datOut.vecXPicardLin = vecXPicardLin;
	end
	%
	if (considerPicard&&considerScl)
		vecTemp = -( matD0 \ (matXF * vecF0) );
		fTemp = vecTemp' * matH0 * vecTemp;
		assert( 0.0 < fTemp );
		vecDeltaPicardScl = vecTemp * (-(vecTemp'*vecG0)/fTemp);
		vecXPicardSclLin  = vecX0 + vecDeltaPicardScl;
		%
		datOut.vecDeltaPicardScl = vecDeltaPicardScl;
		datOut.vecXPicardSclLin = vecXPicardSclLin;
	end
	%
	clear fTemp;
	clear vecTemp;
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
	if (considerPicard)
		%vecXPicardOmega = ...
	end
	%
	if (consdierPicard&&considerScl)
		%vecXPicardSclOmega = ...
	end
	%
	%vecXLevenbergOmega = ...
	%
	if (consdierScl)
		%vecXLevenbergSclOmega = ...
	end
	%
	if (consdierGradCurve)
		%vecXGradCurveOmega = ...
	end
	%
	if (consdierGradCurve&&considerScl)
		%vecXGradCurveSclOmega = ...
	end
	%
	%
	vecXSuggested = vecXNewtonLin;
	%
	%
return;
end

%!test
%!	test_studyPt;
