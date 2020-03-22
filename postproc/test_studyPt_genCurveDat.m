warning("Work-in-progress!");
warning("This should maybe be deprecated in favor of studyPt?");
clear;
commondefs;
thisFile = "test_studyPt_genCurveDat";
tic();
numFigs = 0;
%
sizeX = 9;
sizeF = 7;
sizeK = 5;
%
seedPrm = demoFunc0101_genSeedPrm("moderate");
%randState = mod(round(time),1E6);
%randState = 0;
randState = 677832;
echo__randStat = randState
seedPrm.randState = randState;
seedPrm.sizeX = sizeX;
seedPrm.sizeF = sizeF;
funcPrm = demoFunc0101_genFuncPrm(seedPrm);
%
%
funchF = @(vecXDummy)( demoFunc0101_eval( vecXDummy, funcPrm ) );
funchJ = @(vecXDummy)( demoFunc0101_evalJaco( vecXDummy, funcPrm ) );
%
vecX0 = zeros(seedPrm.sizeX,1);
vecF0 = funchF(vecX0);
matJ = funchJ(vecX0);
matV = eye( sizeX, sizeK );
matW = matJ * matV;
matH = matW' * matW;
vecG = -matW' * vecF0;
vecXSecret = funcPrm.x0;
%
prm = [];
[ curveDat, retCode, datOut ] = studyPt_genCurveDat( ...
  funchF, vecX0, matV, matH, vecG, STEPTYPE__NEWTON, prm );
toc;
%
return;
%
numCurves = size(stepFunchDat.curveDat,2);
assert( issize(stepFunchDat.curveDat,[1,numCurves]) );
for n=1:numCurves
switch (stepFunchDat.curveDat(n).stepType)
case {STEPTYPE__NEWTON}
	stepFunchDat.curveDat(n).col = [ 0.7, 0.0, 0.0 ];
case {STEPTYPE__PICARD}
	stepFunchDat.curveDat(n).col = [ 1.0, 0.0, 1.0 ];
case {STEPTYPE__PICARD_SCALED}
	stepFunchDat.curveDat(n).col = [ 0.5, 0.0, 0.5 ];
case {STEPTYPE__GRADDIR}
	stepFunchDat.curveDat(n).col = [ 0.0, 0.9, 0.0 ];
case {STEPTYPE__GRADDIR_SCALED}
	stepFunchDat.curveDat(n).col = [ 0.0, 0.5, 0.0 ];
case {STEPTYPE__LEVCURVE}
	stepFunchDat.curveDat(n).col = [ 0.9, 0.9, 0.0 ];
case {STEPTYPE__LEVCURVE_SCALED}
	stepFunchDat.curveDat(n).col = [ 0.5, 0.5, 0.0 ];
case {STEPTYPE__GRADCURVE}
	stepFunchDat.curveDat(n).col = [ 0.0, 0.0, 1.0 ];
case {STEPTYPE__GRADCURVE_SCALED}
	stepFunchDat.curveDat(n).col = [ 0.0, 0.0, 0.5 ];
case {STEPTYPE__SECRET}
	stepFunchDat.curveDat(n).col = [ 0.5, 0.5, 0.5 ];
otherwise
	stepFunchDat.curveDat(n).col = [ 0.0, 0.0, 0.0 ];
end
end
%
rvecNuVals = linspace(0.0,1.0,50);
matX0 = repmat( vecX0, size(rvecNuVals) );
matF0 = repmat( vecF0, size(rvecNuVals) );
funchFLin = @(vecXDummy)( matF0 + (matJ*vecXDummy) );
funchOmega = @(vecXDummy)( 0.5*sum(funchF(vecXDummy).^2,1) );
funchOmegaLin = @(vecXDummy)( 0.5*sum(funchFLin(vecXDummy).^2,1) );
%
%
numFigs++; figure(numFigs);
hold off;
for n=1:numCurves
	if (stepFunchDat.curveDat(n).funchYSupportsMultiArg)
		matY = stepFunchDat.curveDat(n).funchYOfNu(rvecNuVals);
	else
		clear matY;
		for m=1:size(rvecNuVals,2)
			matY(:,m) = stepFunchDat.curveDat(n).funchYOfNu(rvecNuVals(m));
		end
	end
	matDelta = matV * matY;
	rvecDeltaNormVals = sqrt(sum( matDelta.^2, 1 ));
	rvecOmegaVals = funchOmega( matX0 + matDelta );
	col = stepFunchDat.curveDat(n).col;
	plot( rvecNuVals, rvecDeltaNormVals, 'o-', ...
	  'color', col, 'markersize', 2*(numCurves+2-n), 'linewidth', 2 );
	hold on;
end
if (0)
	matDelta = (vecXSecret-vecX0)*rvecNuVals;
	rvecDeltaNormVals = sqrt(sum( matDelta.^2, 1 ));
	rvecOmegaVals = funchOmega( matX0 + matDelta );
	plot( rvecNuVals, rvecDeltaNormVals, 'o-', ...
	  'color', [0.0,0.0,0.0], 'markersize', 2, 'linewidth', 2 );
	hold on;
end
grid on;
hold off;
%
numFigs++; figure(numFigs);
hold off;
for n=1:numCurves
	if (stepFunchDat.curveDat(n).funchYSupportsMultiArg)
		matY = stepFunchDat.curveDat(n).funchYOfNu(rvecNuVals);
	else
		clear matY;
		for m=1:size(rvecNuVals,2)
			matY(:,m) = stepFunchDat.curveDat(n).funchYOfNu(rvecNuVals(m));
		end
	end
	matDelta = matV * matY;
	rvecDeltaNormVals = sqrt(sum( matDelta.^2, 1 ));
	rvecOmegaVals = funchOmega( matX0 + matDelta );
	col = stepFunchDat.curveDat(n).col;
	plot( rvecDeltaNormVals, rvecOmegaVals, 'o-', ...
	  'color', col, 'markersize', 2*(numCurves+2-n), 'linewidth', 2 );
	hold on;
end
if (0)
	matDelta = (vecXSecret-vecX0)*rvecNuVals;
	rvecDeltaNormVals = sqrt(sum( matDelta.^2, 1 ));
	rvecOmegaVals = funchOmega( matX0 + matDelta );
	plot( rvecDeltaNormVals, rvecOmegaVals, 'o-', ...
	  'color', [0.0,0.0,0.0], 'markersize', 2, 'linewidth', 2 );
	hold on;
end
grid on;
hold off;
%
numFigs++; figure(numFigs);
hold off;
for n=1:numCurves
	if (stepFunchDat.curveDat(n).funchYSupportsMultiArg)
		matY = stepFunchDat.curveDat(n).funchYOfNu(rvecNuVals);
	else
		clear matY;
		for m=1:size(rvecNuVals,2)
			matY(:,m) = stepFunchDat.curveDat(n).funchYOfNu(rvecNuVals(m));
		end
	end
	matDelta = matV * matY;
	rvecDeltaNormVals = sqrt(sum( matDelta.^2, 1 ));
	rvecOmegaVals = funchOmegaLin( matX0 + matDelta );
	col = stepFunchDat.curveDat(n).col;
	plot( rvecDeltaNormVals, rvecOmegaVals, 'o-', ...
	  'color', col, 'markersize', 2*(numCurves+2-n), 'linewidth', 2 );
	hold on;
end
if (0)
	matDelta = (vecXSecret-vecX0)*rvecNuVals;
	rvecDeltaNormVals = sqrt(sum( matDelta.^2, 1 ));
	rvecOmegaVals = funchOmegaLin( matX0 + matDelta );
	plot( rvecDeltaNormVals, rvecOmegaVals, 'o-', ...
	  'color', [0.0,0.0,0.0], 'markersize', 2, 'linewidth', 2 );
	hold on;
end
grid on;
hold off;
%
return
