clear;
commondefs;
numFigs = 0;
startTime = time();
%
xinit;
funchOmega = @(dummyX)( testfunc2021_funcOmega(dummyX,testFuncPrm) );
funchF = @(x)( testfunc2021_funcF(x,testFuncPrm) );
[ omega0, vecG0, matH0 ] = testfunc2021_funcOmega( vecX0, testFuncPrm );
%
numOptOnF = true;
if (numOptOnF)
	[ vecF0, matJ0 ] = funchF( vecX0, testFuncPrm );
	%echo__omega0 = omega0
	%echo__vecG0 = vecG0
	%echo__matH0 = matH0
	omegaOCQ0 = 0.5*(vecF0'*vecF0)
	vecGOCQ0 = matJ0'*vecF0
	matHOCQ0 = matJ0'*matJ0
	[ matPsiOCQ0, matLambdaOCQ0 ] = eig(matHOCQ0)
	[ lambdaOCQ0AbsMin, nOfOCQ0AbsMin ] = min(abs(diag(matLambdaOCQ0)))
	vecPhiHat = matPsiOCQ0(:,nOfOCQ0AbsMin)
	vecEta = funchF( vecX0 + vecPhiHat ) - ( vecF0 + matJ0*vecPhiHat );
	%
	%funchFOCQ = @(dummyX)( vecF0 + matJ0*(dummyX-vecX0) + vecEta*sumsq(vecPhiHat'*(dummyX-vecX0)) );
	%funchOmegaOCQ = @(dummyX)( sumsq(funchFOCQ(dummyX),1)/2.0 );
	%
	function [ dummyVecF, dummyMatJ ] = funcF_OCQ( dummyX, vecF0, matJ0, vecEta, vecPhiHat, vecX0 )
		dummyVecF = vecF0 + matJ0*(dummyX-vecX0) + vecEta*sumsq(vecPhiHat'*(dummyX-vecX0));
		dummyMatJ = matJ0 + vecEta*vecPhiHat'*(2.0*(vecPhiHat'*(dummyX-vecX0)));
		%
		% The folowing yields csntJ, which is a subset of csntH.
		%dummyVecF = vecF0 + matJ0*(dummyX-vecX0);
		%dummyMatJ = matJ0;
	end
	funchF_OCQ = @(dummyX) funcF_OCQ( dummyX, vecF0, matJ0, vecEta, vecPhiHat, vecX0 );
	%
	function [ dummyOmega, dummyVecG ] = funchOmegaOCQ( dummyX, dummyFuncF )
		[ dummyVecF, dummyMatJ ] = dummyFuncF( dummyX );
		dummyOmega = (dummyVecF'*dummyVecF)/2.0;
		dummyVecG = dummyMatJ'*dummyVecF;
	end
	funchOmega_OCQ = @(dummyX) funchOmegaOCQ( dummyX, funchF_OCQ );
	%
	basicOnF = true;
	if (basicOnF)
		omega0 = omegaOCQ0;
		vecG0 = vecGOCQ0;
		matH0 = matHOCQ0;
	end
end
%
%
%
numCurves = 0;
%
%
if (0)
	% FOR FOCUSED DEV, not curve viz.
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'DEV';
	thisCurvePrm = [];
	thisCurvePrm.matS = [ 5.0, 0.0; 0.0, 1.0 ];
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcFOCQLevCurve( vecX0, vecF0, matJ0, vecPhiHat, vecEta, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
	%
	xpost;
	xviz;
	msg( __FILE__, __LINE__, sprintf("Elapsed time is %0.3fs.", time()-startTime) );
	return
end
%
%
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcGradCurve';
	thisCurvePrm = [];
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcGradCurve( vecX0, funchOmega, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcLevCurve';
	thisCurvePrm = [];
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcLevCurve( vecX0, funchOmega, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcMinfordCurve';
	thisCurvePrm = [];
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcMinfordCurve( vecX0, funchOmega, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
%
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcGradCurve OCQ';
	thisCurvePrm = [];
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcGradCurve( vecX0, funchOmega_OCQ, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcLevCurve OCQ';
	thisCurvePrm = [];
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcLevCurve( vecX0, funchOmega_OCQ, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcMinfordCurve OCQ';
	thisCurvePrm = [];
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcMinfordCurve( vecX0, funchOmega_OCQ, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcFOCQLevCurve';
	thisCurvePrm = [];
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcFOCQLevCurve( vecX0, vecF0, matJ0, vecPhiHat, vecEta, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
if (0)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcBasicLevCurve';
	thisCurvePrm = [];
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcBasicLevCurve( vecX0, omega0, vecG0, matH0, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
if (0)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcBasicGradCurve';
	thisCurvePrm = [];
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcBasicGradCurve( vecX0, omega0, vecG0, matH0, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
xpost;
xviz;
msg( __FILE__, __LINE__, sprintf("Elapsed time is %0.3fs.", time()-startTime) );
