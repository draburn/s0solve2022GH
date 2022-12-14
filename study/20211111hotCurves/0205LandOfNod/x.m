% Started as x3 from 0117lostinparadise.
% This is for "settling" numOpt code from "study" to "numopt";
%  some functionality may be broken.
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
%
doGrad_alytG = false; doLev_dispena = false; doLev_minford = false;
doGrad_cnstH = false; doLev_cnstH = false;
doFOCQ_L = false; doFOCQ_C = false; doFOCQ_R = false;
%
doGrad_alytG = true; doLev_dispena = true; doLev_minford = true;
%doGrad_cnstH = true; doLev_cnstH = true;
doFOCQ_cts = true;
doFOCQ_L = true; doFOCQ_C = true; doFOCQ_R = true;
%
%
numOptOnF = true;
if (numOptOnF)
	[ vecF0, matJ0 ] = funchF( vecX0, testFuncPrm );
	%echo__omega0 = omega0
	%echo__vecG0 = vecG0
	%echo__matH0 = matH0
	omegaOCQ0 = 0.5*(vecF0'*vecF0);
	vecGOCQ0 = matJ0'*vecF0;
	matHOCQ0 = matJ0'*matJ0;
	[ matPsiOCQ0, matLambdaOCQ0 ] = eig(matHOCQ0);
	[ lambdaOCQ0AbsMin, nOfOCQ0AbsMin ] = min(abs(diag(matLambdaOCQ0)));
	vecPhiHat = matPsiOCQ0(:,nOfOCQ0AbsMin);
	%vecEta = ( funchF(vecX0+vecPhiHat) + funchF(vecX0-vecPhiHat) )/2.0 - vecF0
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
doVizFOCQRoots = true;
if (doVizFOCQRoots)
	temp_matPsi = orth( eye(sizeX,sizeX) - (vecPhiHat*(vecPhiHat')), sqrt(eps) );
	temp_vecLambda = matJ0*vecPhiHat;
	temp_matW = matJ0*temp_matPsi;
	prm_vizFOCQRoots = [];
	prm_vizFOCQRoots.figNum = 20;
	vizFOCQRoots( vecF0, temp_vecLambda, vecEta, temp_matW, prm_vizFOCQRoots );
	clear temp_vecLambda;
	clear temp_matPsi;
end
%
%
%
numCurves = 0;
%
%
%
if (doGrad_alytG)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcGradCurve_alytG';
	thisCurvePrm = [];
	curveDat(numCurves).matX = calcGradCurve_alytG( vecX0, funchOmega, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
if (doLev_dispena)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcLevCurve_alytG_dispena';
	thisCurvePrm = [];
	curveDat(numCurves).matX = calcLevCurve_alytG_dispena( vecX0, funchOmega, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
if (doLev_minford)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcLevCurve_alytG_minford';
	thisCurvePrm = [];
	curveDat(numCurves).matX = calcLevCurve_alytG_minford( vecX0, funchOmega, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
%
if (doFOCQ_cts)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcLevCurve_focq cts';
	thisCurvePrm = [];
	thisCurvePrm.curveSelector = 2;
	curveDat(numCurves).matX = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
%
if (doFOCQ_L)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcLevCurve_focq L';
	thisCurvePrm = [];
	thisCurvePrm.curveSelector = -1;
	curveDat(numCurves).matX = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
if (doFOCQ_C)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcLevCurve_focq C';
	thisCurvePrm = [];
	thisCurvePrm.curveSelector = 0;
	curveDat(numCurves).matX = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
if (doFOCQ_R)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcLevCurve_focq R';
	thisCurvePrm = [];
	thisCurvePrm.curveSelector = 1;
	curveDat(numCurves).matX = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhiHat, vecEta, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
%
%
if (doGrad_cnstH)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcGradCurve_cnstH';
	thisCurvePrm = [];
	curveDat(numCurves).matX = calcGradCurve_cnstH( vecX0, omega0, vecG0, matH0, thisCurvePrm );
	thisCurveElapsedTime = time()-thisCurveTime0;
	msg( __FILE__, __LINE__, sprintf( "Calculation of %s took %0.3fs.", thisCurveName, thisCurveElapsedTime ) );
	curveDat(numCurves).elapsedTime = thisCurveElapsedTime;
	curveDat(numCurves).strName = thisCurveName;
	curveDat(numCurves).prm = thisCurvePrm;
end
if (doLev_cnstH)
	numCurves++;
	thisCurveTime0 = time();
	thisCurveName = 'calcLevCurve_cnstH';
	thisCurvePrm = [];
	curveDat(numCurves).matX = calcLevCurve_cnstH( vecX0, omega0, vecG0, matH0, thisCurvePrm );
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
