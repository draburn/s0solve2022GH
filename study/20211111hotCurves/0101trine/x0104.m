clear;
thisFile = "x0104";
commondefs;
numFigs = 0;
startTime = time();
%
xinit; thisFile = "x0104";
funchF = @(x)( testFunc_eval(x,testFuncPrm) );
funchOmega = @(x)( 0.5*sum(testFunc_eval(x,testFuncPrm).^2) );
funchG = @(x)( testFunc_evalG(x,testFuncPrm) );
funchGOmega = @(x)( testFunc_evalGOmega(x,testFuncPrm) );%
%
%
numCurves = 0;
%
%
if (0)
	% FOR FOCUSED DEV, not curve viz.
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = '<<DEV>>';
	%matS = [];
	matS = [ 5.0, 0.0; 0.0, 1.0 ];
	thisCurvePrm = [];
	msg( thisFile, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcMinfordCurve( funchOmega, funchG, vecX0, matS, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( thisFile, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
%xpost; thisFile = "x1";
%xviz; thisFile = "x1";
msg( thisFile, __LINE__, sprintf("Elapsed time is %0.3fs.", time()-startTime) );
thisFile = "RETURN FROM x0104";
return
end
%
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = "calcGradescentCurve_lsode";
	thisCurvePrm = [];
	msg( thisFile, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcGradescentCurve_lsode( funchOmega, funchG, vecX0, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( thisFile, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcLevCurve (fminunc?)';
	thisCurvePrm = [];
	msg( thisFile, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcLevCurve( funchOmega, funchG, vecX0, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( thisFile, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcMinfordCurve_adaptiveDeltaR';
	matS = [];
	thisCurvePrm = [];
	msg( thisFile, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcMinfordCurve_adaptiveDeltaR( funchOmega, funchG, vecX0, matS, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( thisFile, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcMinfordCurve';
	matS = [];
	thisCurvePrm = [];
	msg( thisFile, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcMinfordCurve( funchOmega, funchG, vecX0, matS, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( thisFile, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcMinfordCurve (scaled [5,0;0,1])';
	matS = [ 5.0, 0.0; 0.0, 1.0 ];
	thisCurvePrm = [];
	msg( thisFile, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcMinfordCurve( funchOmega, funchG, vecX0, matS, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( thisFile, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
	% This is not a comment, it is to force a desired git diff.
end
%
%
if (1)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcMinfordCurve (scaled [1,0;0,5])';
	matS = [ 1.0, 0.0; 0.0, 5.0 ];
	thisCurvePrm = [];
	msg( thisFile, __LINE__, sprintf( "Calculating %s...", thisCurveName ) );
	curveDat(numCurves).matX = calcMinfordCurve( funchOmega, funchG, vecX0, matS, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( thisFile, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
	% This is not a comment, it is to force a desired git diff.
end
%
xpost; thisFile = "x0104";
xviz; thisFile = "x0104";
msg( thisFile, __LINE__, sprintf("Elapsed time is %0.3fs.", time()-startTime) );
thisFile = "RETURN FROM x0104";
return;
