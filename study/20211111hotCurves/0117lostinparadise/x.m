clear;
commondefs;
numFigs = 0;
startTime = time();
%
xinit;
funchOmega = @(dummyX)( testfunc2021_funcOmega(dummyX,testFuncPrm) );
funchF = @(x)( testfunc2021_funcF(x,testFuncPrm) );
[ omega0, vecG0, matH0 ] = testfunc2021_funcOmega( vecX0, testFuncPrm );
matS_shared = [ 1.0, 0.0; 0.0, 3.0 ];
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
	curveDat(numCurves).matX = calcBasicGradCurve( vecX0, omega0, vecG0, matH0, thisCurvePrm );
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
	thisCurveName = 'calcGradCurve Scaled';
	thisCurvePrm = [];
	thisCurvePrm.matS = matS_shared;
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
	thisCurveName = 'calcLevCurve Scaled';
	thisCurvePrm = [];
	thisCurvePrm.matS = matS_shared;
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
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcMinfordCurve Scaled';
	thisCurvePrm = [];
	thisCurvePrm.matS = matS_shared;
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcMinfordCurve( vecX0, funchOmega, thisCurvePrm );
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
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcBasicLevCurve Scaled';
	thisCurvePrm = [];
	thisCurvePrm.matS = matS_shared;
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcBasicLevCurve( vecX0, omega0, vecG0, matH0, thisCurvePrm );
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
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcBasicGradCurve Scaled';
	thisCurvePrm = [];
	thisCurvePrm.matS = matS_shared;
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
