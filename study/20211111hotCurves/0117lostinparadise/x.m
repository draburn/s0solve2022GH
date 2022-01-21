clear;
commondefs;
numFigs = 0;
startTime = time();
%
xinit;
funchOmega = @(dummyX)( testfunc2021_funcOmega(dummyX,testFuncPrm) );
funchF = @(x)( testfunc2021_funcF(x,testFuncPrm) );
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
	%matS = [];
	matS = [ 5.0, 0.0; 0.0, 1.0 ];
	thisCurvePrm = [];
	msg( __FILE__, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcMinfordCurve( vecX0, funchOmega, thisCurvePrm );
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
	%matS = [];
	matS = [ 5.0, 0.0; 0.0, 1.0 ];
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
	matS = [];
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
	matS = [];
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
xpost;
xviz;
msg( __FILE__, __LINE__, sprintf("Elapsed time is %0.3fs.", time()-startTime) );
