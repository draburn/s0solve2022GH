numFigs = 0;
%
numPts = 101;
pPts = linspace( 0.0, 1.0, numPts );
vecDeltaPts = zeros(sizeX,numPts);
vecFModelPts = zeros(sizeF,numPts);
vecFPts = zeros(sizeF,numPts);
omegaHModelPts = zeros(1,numPts);
for n=1:numPts
	vecDelta = funchDeltaOfP( pPts(n) );
	vecDeltaPts(:,n) = vecDelta;
	vecFModelPts(:,n) = funchFModelOfDelta( vecDelta );
	vecFPts(:,n) = funchF( vecX + vecDelta );
	omegaHModelPts(n) = funchOmegaHessModelOfDelta( vecDelta );
endfor
dPts = sqrt(sumsq(vecDeltaPts,1));
omegaFModelPts = sumsq(vecFModelPts,1)/2.0;
omegaPts = sumsq(vecFPts,1)/2.0;
%
numFigs++; figure(numFigs);
plot( ...
  dPts, omegaHModelPts, 'o-', ...
  dPts, omegaFModelPts, 'x-', ...
  dPts, omegaPts, 's-' );
grid on;
legend( ...
  "omegaHModel", ...
  "omegaFModel", ...
  "omega", ...
  "location", "northeast" );
return;
