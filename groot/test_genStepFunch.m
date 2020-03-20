clear;
commondefs;
thisFile = "test_genStepFunch";
tic();
numFigs = 0;
%
sizeX = 10;
sizeF = 10;
sizeK = 10;
%
seedPrm = demoFunc0101_genSeedPrm("moderate");
%randState = mod(round(time),1E6)
randState = 0;
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
%
prm = [];
[ retCode, studyPtDat ] = genStepFunch( funchF, vecX0, matW, matV, prm );
%
rvecNuVals = linspace(0.0,1.0,100);
matX0 = repmat( vecX0, size(rvecNuVals) );
matF0 = repmat( vecF0, size(rvecNuVals) );
funchFLin = @(vecXDummy)( matF0 + (matJ*vecXDummy) );
%funchOmega = @(vecXDummy)( 0.5*sum(funchF(vecXDummy).^2,1) );
funchOmega = @(vecXDummy)( 0.5*sum(funchFLin(vecXDummy).^2,1) );
%
numCurves = size(studyPtDat.curveDat,2)
%
numFigs++; figure(numFigs);
hold off;
for n=1:numCurves
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
	rvecOmegaVals = funchOmega( matX0 + matDelta );
	plot( rvecDeltaNormVals, rvecOmegaVals, 'o-' );
	hold on;
	grid on;
end
hold off;
%
return
