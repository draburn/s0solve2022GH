clear;
thisFile = "x1";
commondefs;
numFigs = 0;
startTime = time();
%
xinit; thisFile = "x1";
funchF = @(x)( testFunc_eval(x,testFuncPrm) );
funchOmegaG = @(x)( testFunc_evalOmegaG(x,testFuncPrm) );
funchGOmega = @(x)( testFunc_evalGOmega(x,testFuncPrm) );
funchGForLSODE = @(x,t)( testFunc_evalG(x,testFuncPrm) );
%
funchGradDescent = @(x)( -testFunc_evalG(x,testFuncPrm) );
funchG = @(x)( testFunc_evalG(x,testFuncPrm) );
%
%
numCurves = 0;
%
%
time0 = time();
numCurves++;
curveDat(numCurves).matX = NEO_calcGradCurve1128( funchGForLSODE, vecX0 );
curveDat(numCurves).strName = 'NEO\_calcGradCurve1128';
msg( thisFile, __LINE__, sprintf( ...
  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
%
%
if (1)
	time0 = time();
	numCurves++;
	curveDat(numCurves).matX = calcGradDescentCurve_alpha( funchGOmega, vecX0 );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_alpha';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );

end
%
%
if (0)
	time0 = time();
	numCurves++;
	curveDat(numCurves).matX = OLD_calcHOTLevCurve1112( testFuncPrm, vecX0 );
	curveDat(numCurves).strName = 'OLD\_calcHOTLevCurve1112';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
%
%
if (0)
	time0 = time();
	numCurves++;
	curveDat(numCurves).matX = NEO_calcLevCurve_lsodeBarrier( funchGForLSODE, vecX0 );
	curveDat(numCurves).strName = 'NEO\_calcLevCurve\_lsodeBarrier';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
%
%
if (0)
	time0 = time();
	numCurves++;
	curveDat(numCurves).matX = NEO_calcLevCurve_lsodePlusCost( funchGForLSODE, vecX0 );
	curveDat(numCurves).strName = 'NEO\_calcLevCurve\_lsodePlusCost';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
%
%
if (0)
	time0 = time();
	numCurves++;
	curveDat(numCurves).matX = NEO_calcLevCurve_lsodeOnSurface( funchGForLSODE, vecX0 );
	curveDat(numCurves).strName = 'NEO\_calcLevCurve\_lsodeOnSurface';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
%
%
if (0)
	time0 = time();
	numCurves++;
	curveDat(numCurves).matX = NEO_calcLevCurve_general( funchGForLSODE, vecX0 );
	curveDat(numCurves).strName = 'NEO\_calcLevCurve\_general';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
%
%
if (0)
	time0 = time();
	numCurves++;
	curveDat(numCurves).matX = calcGradDescentCurve( funchGradDescent, vecX0 );
	curveDat(numCurves).strName = 'calcGradDescentCurve';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
%
%
%
xpost; thisFile = "x1";
xviz; thisFile = "x1";
msg( thisFile, __LINE__, sprintf("Elapsed time is %0.3fs.", time()-startTime) );
thisFile = "RETURN FROM x1";
return;
