clear;
commondefs;
thisFile = "test_studyPt";
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
%randState = 677832;
randState = 686206;
echo__randState = randState
randn("seed",randState);
rand("seed",randState);
seedPrm.randState = randState;
seedPrm.sizeX = sizeX;
seedPrm.sizeF = sizeF;
funcPrm = demoFunc0101_genFuncPrm(seedPrm);
randn("seed",randState);
rand("seed",randState);
%
%
funchF = @(vecXDummy)( demoFunc0101_eval( vecXDummy, funcPrm ) );
funchJ = @(vecXDummy)( demoFunc0101_evalJaco( vecXDummy, funcPrm ) );
%
vecX0 = zeros(seedPrm.sizeX,1);
vecF0 = funchF(vecX0);
matJ = funchJ(vecX0);
%matV = eye( sizeX, sizeK );
matU = randn(sizeX,sizeK);
matV = myorth(matU);
matW = matJ * matV;
vecXSecret = funcPrm.x0;
%
prm = [];
[ retCode, studyPtDat ] = studyPt( ...
  funchF, vecX0, matW, matV, vecXSecret );
%studyPtDat.curveDat = studyPtDat.curveDat([8]);
%
numCurves = size(studyPtDat.curveDat,2);
assert( issize(studyPtDat.curveDat,[1,numCurves]) );
for n=1:numCurves
switch (studyPtDat.curveDat(n).stepType)
case {STEPTYPE__NEWTON}
	studyPtDat.curveDat(n).col = [ 0.7, 0.0, 0.0 ];
case {STEPTYPE__PICARD}
	studyPtDat.curveDat(n).col = [ 1.0, 0.0, 1.0 ];
case {STEPTYPE__PICARD_SCALED}
	studyPtDat.curveDat(n).col = [ 0.5, 0.0, 0.5 ];
case {STEPTYPE__GRADDIR}
	studyPtDat.curveDat(n).col = [ 0.0, 0.9, 0.0 ];
case {STEPTYPE__GRADDIR_SCALED}
	studyPtDat.curveDat(n).col = [ 0.0, 0.5, 0.0 ];
case {STEPTYPE__LEVCURVE}
	studyPtDat.curveDat(n).col = [ 0.9, 0.9, 0.0 ];
case {STEPTYPE__LEVCURVE_SCALED}
	studyPtDat.curveDat(n).col = [ 0.5, 0.5, 0.0 ];
case {STEPTYPE__GRADCURVE}
	studyPtDat.curveDat(n).col = [ 0.0, 0.0, 1.0 ];
case {STEPTYPE__GRADCURVE_SCALED}
	studyPtDat.curveDat(n).col = [ 0.0, 0.0, 0.5 ];
case {STEPTYPE__SECRET}
	studyPtDat.curveDat(n).col = [ 0.5, 0.5, 0.5 ];
otherwise
	studyPtDat.curveDat(n).col = [ 0.0, 0.0, 0.0 ];
end
end
%
numNuVals = 100;
funchFLin = @(vecXDummy)( repmat(vecF0,[1,size(vecXDummy,2)]) + (matJ*vecXDummy) );
funchOmega = @(vecXDummy)( 0.5*sum(funchF(vecXDummy).^2,1) );
funchOmegaLin = @(vecXDummy)( 0.5*sum(funchFLin(vecXDummy).^2,1) );
%
%
numFigs++; figure(numFigs);
hold off;
for n=1:numCurves
	rvecNuVals = studyPtDat.curveDat(n).rvecNuVals;
	if (studyPtDat.curveDat(n).funchYSupportsMultiArg)
		matY = studyPtDat.curveDat(n).funchYOfNu(rvecNuVals);
	else
		clear matY;
		for m=1:size(rvecNuVals,2)
			matY(:,m) = studyPtDat.curveDat(n).funchYOfNu(rvecNuVals(m));
		end
	end
	matDelta = matV * matY;
	rvecDeltaNormVals = sqrt(sum( matDelta.^2, 1 ));
	rvecOmegaVals = funchOmega( repmat(vecX0,[1,size(matDelta,2)]) + matDelta );
	col = studyPtDat.curveDat(n).col;
	plot( rvecNuVals, rvecDeltaNormVals, 'o-', ...
	  'color', col, 'markersize', 2*(numCurves+2-n), 'linewidth', 2 );
	hold on;
end
if (0)
	rvecNuVals = linspace(0.0,1.0,100);
	matDelta = (vecXSecret-vecX0)*rvecNuVals;
	rvecDeltaNormVals = sqrt(sum( matDelta.^2, 1 ));
	rvecOmegaVals = funchOmega( repmat(vecX0,[1,size(matDelta,2)]) + matDelta );
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
	rvecNuVals = studyPtDat.curveDat(n).rvecNuVals;
	if (studyPtDat.curveDat(n).funchYSupportsMultiArg)
		matY = studyPtDat.curveDat(n).funchYOfNu(rvecNuVals);
	else
		clear matY;
		for m=1:size(rvecNuVals,2)
			matY(:,m) = studyPtDat.curveDat(n).funchYOfNu(rvecNuVals(m));
		end
	end
	matDelta = matV * matY;
	rvecDeltaNormVals = sqrt(sum( matDelta.^2, 1 ));
	rvecOmegaVals = funchOmega( repmat(vecX0,[1,size(matDelta,2)]) + matDelta );
	col = studyPtDat.curveDat(n).col;
	plot( rvecDeltaNormVals, rvecOmegaVals, 'o-', ...
	  'color', col, 'markersize', 2*(numCurves+2-n), 'linewidth', 2 );
	hold on;
end
if (0)
	rvecNuVals = linspace(0.0,1.0,100);
	matDelta = (vecXSecret-vecX0)*rvecNuVals;
	rvecDeltaNormVals = sqrt(sum( matDelta.^2, 1 ));
	rvecOmegaVals = funchOmega( repmat(vecX0,[1,size(matDelta,2)]) + matDelta );
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
	rvecNuVals = studyPtDat.curveDat(n).rvecNuVals;
	if (studyPtDat.curveDat(n).funchYSupportsMultiArg)
		matY = studyPtDat.curveDat(n).funchYOfNu(rvecNuVals);
	else
		clear matY;
		for m=1:size(rvecNuVals,2)
			matY(:,m) = studyPtDat.curveDat(n).funchYOfNu(rvecNuVals(m));
		end
	end
	matDelta = matV * matY;
	rvecDeltaNormVals = sqrt(sum( matDelta.^2, 1 ));
	rvecOmegaVals = funchOmegaLin( repmat(vecX0,[1,size(matDelta,2)]) + matDelta );
	col = studyPtDat.curveDat(n).col;
	plot( rvecDeltaNormVals, rvecOmegaVals, 'o-', ...
	  'color', col, 'markersize', 2*(numCurves+2-n), 'linewidth', 2 );
	hold on;
end
if (0)
	rvecNuVals = linspace(0.0,1.0,100);
	matDelta = (vecXSecret-vecX0)*rvecNuVals;
	rvecDeltaNormVals = sqrt(sum( matDelta.^2, 1 ));
	rvecOmegaVals = funchOmegaLin( repmat(vecX0,[1,size(matDelta,2)]) + matDelta );
	plot( rvecDeltaNormVals, rvecOmegaVals, 'o-', ...
	  'color', [0.0,0.0,0.0], 'markersize', 2, 'linewidth', 2 );
	hold on;
end
grid on;
hold off;
%
return
