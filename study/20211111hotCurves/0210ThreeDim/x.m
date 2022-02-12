clear;
commondefs;
numFigs = 0;
startTime = time();
%
%
xinit;
%
doGrad_alytG = true; doLev_dispena = true; doLev_minford = true;
doFOCQ_minXi0_jtj = true;
doFOCQ_minXi0_fullish = true;
%doFOCQ_L0jtj = true; doFOCQ_C0jtj = true; doFOCQ_R0jtj = true;
doGrad_cnstH_jtj = true; doLev_cnstH_jtj = true; doGradSeg_cnstH_jtj = true;
doGrad_cnstH_full = true; doLev_cnstH_full = true; doGradSeg_cnstH_full = true;
%
msg( __FILE__, __LINE__, "Calculating curves..." );
numCurves = 0;
%
%
if ( doGrad_alytG )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcGradCurve_alytG';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).vecXVals = calcGradCurve_alytG( vecX0, funchOmega, curveDat(numCurves).prm );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
if ( doLev_dispena )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcLevCurve_alytG_dispena';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).vecXVals = calcLevCurve_alytG_dispena( vecX0, funchOmega, curveDat(numCurves).prm );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
if ( doLev_minford )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcLevCurve_alytG_minford';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).vecXVals = calcLevCurve_alytG_minford( vecX0, funchOmega, curveDat(numCurves).prm );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
%
%
%
if ( doFOCQ_minXi0_jtj )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcLevCurve_focq minXi0 jtj';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).curveSelector = -2;
	curveDat(numCurves).vecXVals = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhi0_jtj, vecEta0_jtj, curveDat(numCurves).prm  );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
if ( doFOCQ_minXi0_fullish )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcLevCurve_focq minXi0 full(ish)';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).curveSelector = -2;
	curveDat(numCurves).vecXVals = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhi0_full, vecEta0_full, curveDat(numCurves).prm  );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
%
%
if ( doFOCQ_L0jtj )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcLevCurve_focq L0jtj';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).curveSelector = -1;
	curveDat(numCurves).vecXVals = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhi0_jtj, vecEta0_jtj, curveDat(numCurves).prm  );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
if ( doFOCQ_C0jtj )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcLevCurve_focq C0jtj';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).curveSelector = 0;
	curveDat(numCurves).vecXVals = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhi0_jtj, vecEta0_jtj, curveDat(numCurves).prm  );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
if ( doFOCQ_R0jtj )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcLevCurve_focq R0jtj';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).curveSelector = 1;
	curveDat(numCurves).vecXVals = calcLevCurve_focq( vecX0, vecF0, matJ0, vecPhi0_jtj, vecEta0_jtj, curveDat(numCurves).prm  );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
%
%
if ( doGrad_cnstH_jtj )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcGradCurve_cnstH jtj';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).vecXVals = calcGradCurve_cnstH( vecX0, omega0, vecG0, matH0_jtj, curveDat(numCurves).prm );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
if ( doLev_cnstH_jtj )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcLevCurve_cnstH jtj';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).vecXVals = calcLevCurve_cnstH( vecX0, omega0, vecG0, matH0_jtj, curveDat(numCurves).prm );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
if ( doGradSeg_cnstH_jtj )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcGradSeg_cnstH jtj';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).vecXVals = calcGradSeg_cnstH( vecX0, omega0, vecG0, matH0_jtj, curveDat(numCurves).prm );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
%
%
if ( doGrad_cnstH_full )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcGradCurve_cnstH full';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).vecXVals = calcGradCurve_cnstH( vecX0, omega0, vecG0, matH0_full, curveDat(numCurves).prm );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
if ( doLev_cnstH_full )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcLevCurve_cnstH full';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).vecXVals = calcLevCurve_cnstH( vecX0, omega0, vecG0, matH0_full, curveDat(numCurves).prm );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
if ( doGradSeg_cnstH_full )
	numCurves++;
	curveDat(numCurves).startTime = time();
	curveDat(numCurves).strName = 'calcGradSeg_cnstH full';
	curveDat(numCurves).prm = [];
	curveDat(numCurves).vecXVals = calcGradSeg_cnstH( vecX0, omega0, vecG0, matH0_full, curveDat(numCurves).prm );
	curveDat(numCurves).numPts = size( curveDat(numCurves).vecXVals, 2 );
	curveDat(numCurves).elapsedTime = time() - curveDat(numCurves).startTime;
	msg( __FILE__, __LINE__, sprintf( "  %s: %d pts in %0.3fs.", ...
	  curveDat(numCurves).strName, curveDat(numCurves).numPts, curveDat(numCurves).elapsedTime ) );
end
%
%
%
xpost;
xviz;
msg( __FILE__, __LINE__, sprintf("Elapsed time is %0.3fs.", time()-startTime) );
