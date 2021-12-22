clear;
thisFile = "x1";
commondefs;
numFigs = 0;
startTime = time();
%
xinit; thisFile = "x1";
funchF = @(x)( testFunc_eval(x,testFuncPrm) );
funchOmega = @(x)( 0.5*sum(testFunc_eval(x,testFuncPrm).^2) );
funchG = @(x)( testFunc_evalG(x,testFuncPrm) );
funchGOmega = @(x)( testFunc_evalGOmega(x,testFuncPrm) );
%%%funchXDot_lsode = @(x,t)( -testFunc_evalG(x',testFuncPrm)' );
%%%funchGradDescent = @(x)( -testFunc_evalG(x,testFuncPrm) );

if (1)
	msg( thisFile, __LINE__, "~~~ Retro HACK 2021.12.03.0100! ~~~~" );
	funchGForLSODE = @(x,t)( testFunc_evalG(x,testFuncPrm) );
	addpath( "study/20211111hotCurves/1128refine" );
	addpath( "study/20211111hotCurves/1201ballin" ); % Make me first!
end


%
%
%
numCurves = 0;
%
%
if (0)
	time0 = time();
	numCurves++;
	myprm = [];
	myprm.stepOrder = 1;
	curveDat(numCurves).matX = calcGradDescentCurve_blind( funchOmega, funchG, vecX0, myprm );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_blind (Euler)';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
if (0)
	time0 = time();
	numCurves++;
	myprm = [];
	myprm.stepOrder = -1;
	curveDat(numCurves).matX = calcGradDescentCurve_blind( funchOmega, funchG, vecX0, myprm );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_blind (Adams-Bashforth)';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
%
%
if (0)
	time0 = time();
	numCurves++;
	myprm = [];
	curveDat(numCurves).matX = calcGradDescentCurve_crudeMomentum( funchOmega, funchG, vecX0, myprm );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_crudeMomentum';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
%
%
if (0)
	time0 = time();
	numCurves++;
	myprm = [];
	myprm.stepOrder = 1;
	curveDat(numCurves).matX = calcGradDescentCurve_crude( funchOmega, funchG, vecX0, myprm );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_crude (1)';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
if (0)
	time0 = time();
	numCurves++;
	myprm = [];
	myprm.stepOrder = 2;
	curveDat(numCurves).matX = calcGradDescentCurve_crude( funchOmega, funchG, vecX0, myprm );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_crude (2)';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
if (0)
	time0 = time();
	numCurves++;
	myprm = [];
	myprm.stepOrder = 4;
	curveDat(numCurves).matX = calcGradDescentCurve_crude( funchOmega, funchG, vecX0, myprm );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_crude (4)';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
%
%
if (1)
	time0 = time();
	numCurves++;
	cgdcbPrm = [];
	cgdcbPrm.doLev = false;
	cgdcbPrm.stepSize = 0.5; % Large, to emphasize inaccuracy.
	curveDat(numCurves).matX = calcGradDescentCurve_beta( funchOmega, funchG, vecX0, cgdcbPrm );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_beta (grad coarse)';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
if (0)
	time0 = time();
	numCurves++;
	cgdcbPrm = [];
	cgdcbPrm.doLev = false;
	cgdcbPrm.stepSize = 0.01;
	cgdcbPrm.iterLimit = 1000;
	curveDat(numCurves).matX = calcGradDescentCurve_beta( funchOmega, funchG, vecX0, cgdcbPrm );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_beta (grad fine)';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
if (1)
	time0 = time();
	numCurves++;
	cgdcbPrm = [];
	cgdcbPrm.doLev = true;
	curveDat(numCurves).matX = calcGradDescentCurve_beta( funchOmega, funchG, vecX0, cgdcbPrm );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_beta (Lev)';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
%
%
if (1)
	time0 = time();
	numCurves++;
	curveDat(numCurves).matX = calcGradCurve_lsode( funchG, vecX0 );
	curveDat(numCurves).strName = 'calcGradCurve\_lsode';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
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
if (0)
	time0 = time();
	numCurves++;
	[ foo, curveDat(numCurves).matX ] = calcGradDescentCurve_alpha( funchGOmega, vecX0 );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_alpha (every geval)';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
%
%
%
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
	% Sometimes fails?
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
if (0)
	time0 = time();
	numCurves++;
	myprm = [];
	curveDat(numCurves).matX = calcGradDescentCurve_eulerMin( funchOmega, funchG, vecX0, myprm );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_eulerMin';
	msg( thisFile, __LINE__, sprintf( ...
	  "Calculation of curve \"%s\" took %0.3fs.", curveDat(numCurves).strName, time()-time0 ) );
end
%
%
if (1)
	time0 = time();
	numCurves++;
	myprm = [];
	curveDat(numCurves).matX = calcGradCurve_chunky( funchG, vecX0, myprm );
	curveDat(numCurves).strName = 'calcGradDescentCurve\_chunky';
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
