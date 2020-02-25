function [ curveDat, retCode, datOut ] = getLLMCurves( vecF, matW, numPts, prm=[], datIn=[] )
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% BASIC INIT.
	%
	commondefs;
	getLLMCurves_setCnsts;
	thisFile = "getCurves";
	retCode = RETCODE__NOT_SET;
	startTime = time();
	%
	% Verbosity.
	verbLev = mygetfield( prm, "verbLev", VERBLEV__PROGRESS );
	reportInterval = mygetfield( prm, "reportInterval", 0.0 );
	assert( isrealscalar(verbLev) );
	assert( isrealscalar(reportInterval) );
	assert( 0.0 <= reportInterval );
	reportTimePrev = startTime - 0.1;
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% VALIDATE INPUT
	%
	%
	sizeF = size(vecF,1);
	assert( 1 <= sizeF );
	assert( isrealarray(vecF,[sizeF,1]) );
	sizeK = size(matW,2);
	assert( 1 <= sizeK );
	assert( isrealarray(matW,[sizeF,sizeK]) );
	assert( 2 <= numPts );
	%
	curveTypes_default = [ ...
	  GETCURVES_CURVETYPE__NEWTON, ...
	  GETCURVES_CURVETYPE__PICARD, ...
	  GETCURVES_CURVETYPE__PICARD_SCALED, ...
	  GETCURVES_CURVETYPE__GRADDIR, ...
	  GETCURVES_CURVETYPE__GRADDIR_SCALED, ...
	  GETCURVES_CURVETYPE__LEVCURVE, ...
	  GETCURVES_CURVETYPE__LEVCURVE_SCALED, ...
	  GETCURVES_CURVETYPE__GRADCURVE, ...
	  GETCURVES_CURVETYPE__GRADCURVE_SCALED ];
	curveTypes = mygetfield( prm, "curveTypes", curveTypes_default );
	numCurves = max(size(curveTypes));
	assert( 1 <= numCurves );
	%
	matV = mygetfield( prm, "matV", [] );
	if (isempty(matV))
		matV = eye(sizeF,sizeK);
	end
	sizeX = size(matV,1);
	assert( 1 <= sizeX );
	assert( isrealarray(matV,[sizeX,sizeK]) );
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% DO WORK.
	%
	matI = eye(sizeK,sizeK);
	vecG = matW' * vecF;
	matH = matW' * matW;
	vecD = diag(matH);
	assert( min(vecD) > 0.0 );
	matD = diag(diag(matH));
	sVals = (0:numPts-1)/(numPts-1.0);
	%
	for curveIndex = 1 : numCurves
		switch (curveTypes(curveIndex))
		case {GETCURVES_CURVETYPE__NEWTON}
			vecTemp = -matH\vecG;
			denom = vecTemp'*matH*vecTemp;
			assert( denom > 0.0 );
			s0 = -(vecTemp'*vecG)/denom;
			assert( abs(s0-1.0) < (sizeK^3)*10.0*(eps^0.75) );
			vecTemp *= s0;
			curveDat(curveIndex).matY = vecTemp * sVals;
			curveDat(curveIndex).strType = "Newton";
			clear denom;
			clear s0;
			clear vecTemp;
		case {GETCURVES_CURVETYPE__PICARD}
			if (sizeF~=sizeK)
				msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
				  "Is prm.matV set properly?" ) );
			end
			vecTemp = -matV' * vecF;
			denom = vecTemp'*matH*vecTemp;
			assert( denom > 0.0 );
			s0 = -(vecTemp'*vecG)/denom;
			vecTemp *= s0;
			curveDat(curveIndex).matY = vecTemp * sVals;
			curveDat(curveIndex).strType = "Picard";
			clear denom;
			clear s0;
			clear vecTemp;
		case {GETCURVES_CURVETYPE__PICARD_SCALED}
			if (sizeF~=sizeK)
				msg_warn( verbLev, thisFile, __LINE__, sprintf( ...
				  "Is prm.matV set properly?" ) );
			end
			vecTemp = -matD\(matV'*vecF);
			denom = vecTemp'*matH*vecTemp;
			assert( denom > 0.0 );
			s0 = -(vecTemp'*vecG)/denom;
			vecTemp *= s0;
			curveDat(curveIndex).matY = vecTemp * sVals;
			curveDat(curveIndex).strType = "PicardScl";
			clear denom;
			clear s0;
			clear vecTemp;
		case {GETCURVES_CURVETYPE__GRADDIR}
			vecTemp = -vecG;
			denom = vecTemp'*matH*vecTemp;
			assert( denom > 0.0 );
			s0 = -(vecTemp'*vecG)/denom;
			vecTemp *= s0;
			curveDat(curveIndex).matY = vecTemp * sVals;
			curveDat(curveIndex).strType = "GradDir";
			clear denom;
			clear s0;
			clear vecTemp;
		case {GETCURVES_CURVETYPE__GRADDIR_SCALED}
			vecTemp = -matD\vecG;
			denom = vecTemp'*matH*vecTemp;
			assert( denom > 0.0 );
			s0 = -(vecTemp'*vecG)/denom;
			vecTemp *= s0;
			curveDat(curveIndex).matY = vecTemp * sVals;
			curveDat(curveIndex).strType = "GradDirScl";
			clear denom;
			clear s0;
			clear vecTemp;
		case {GETCURVES_CURVETYPE__LEVCURVE}
			matL = max(vecD) * matI;
			matA = matH-matL;
			funchY = @(nu)( -nu*((matL+(nu*matA))\vecG) );
			funchF = @(nu)( sqrt(sum((funchY(nu)).^2)) );
			nuVals = flinspace( 0.0, 1.0, numPts, funchF );
			for n=1:max(size(nuVals))
				curveDat(curveIndex).matY(:,n) = funchY(nuVals(n));
			end
			curveDat(curveIndex).strType = "Leveneberg";
			clear nuVals;
			clear funchF;
			clear funchY;
			clear matA;
			clear matL;
		otherwise
			error(sprintf( "Value of curveTypes(%g) is invalid (%g).", ...
			  curveIndex, curveTypes(curveIndex) ));
		end
	end
	%
	retCode = RETCODE__SUCCESS;
	datOut = [];
return;
end

%!test
%!	commondefs;
%!	getLLMCurves_setCnsts;
%!	thisFile = "test getLLMCurvs";
%!	vecF = [1;2]
%!	matW = [1,2;3,5]
%!	numPts = 20;
%!	curveTypes = [ ...
%!	  GETCURVES_CURVETYPE__NEWTON, ...
%!	  GETCURVES_CURVETYPE__PICARD, ...
%!	  GETCURVES_CURVETYPE__PICARD_SCALED, ...
%!	  GETCURVES_CURVETYPE__GRADDIR, ...
%!	  GETCURVES_CURVETYPE__GRADDIR_SCALED ];
%!	prm.curveTypes = curveTypes;
%!	[ curveDat, retCode, datOut ] = getLLMCurves( vecF, matW, numPts, prm );
