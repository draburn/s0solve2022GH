numFigs = 0;
%
numPts = 1001;
pPts = linspace( 0.0, 1.0, numPts );
pPts = 1.0 - ( 1.0 - pPts.^2).^2;
vecDeltaPts = zeros(sizeX,numPts);
vecFModelPts = zeros(sizeF,numPts);
omegaHModelPts = zeros(1,numPts);
for n=1:numPts
	vecDelta = funchDeltaOfP( pPts(n) );
	vecDeltaPts(:,n) = vecDelta;
	vecFModelPts(:,n) = funchFModelOfDelta( vecDelta );
	%%% fevalCount++???
	omegaHModelPts(n) = funchOmegaHessModelOfDelta( vecDelta );
endfor
dPts = sqrt(sumsq(vecDeltaPts,1));
omegaFModelPts = sumsq(vecFModelPts,1)/2.0;
%
%
numPts2 = 51;
pPts2 = linspace( 0.0, 1.0, numPts2 );
pPts2 = 1.0 - ( 1.0 - pPts2.^2).^2;
vecFPts2 = zeros(sizeF,numPts2);
for n=1:numPts2
	vecDelta = funchDeltaOfP( pPts2(n) );
	vecDeltaPts2(:,n) = vecDelta;
	vecFPts2(:,n) = funchF( vecX + vecDelta );
endfor
dPts2 = sqrt(sumsq(vecDeltaPts2,1));
omegaPts2 = sumsq(vecFPts2,1)/2.0;
%
numFigs++; figure(numFigs);
plot( ...
  dPts, omegaHModelPts, 's-', ...
  dPts, omegaFModelPts, 'o-', ...
  dPts2, omegaPts2, 'x-', 'linewidth', 3, 'markersize', 15 );
grid on;
legend( ...
  "omegaHModel", ...
  "omegaFModel", ...
  "omega", ...
  "location", "northeast" );
%
