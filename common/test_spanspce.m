clear;
tic();
sizeX = 5;
sizeK = 2;
randnState = mod(round(time),1E6)
%randnState = 0;
echo_randnState = randnState
randn( "state", randnState );
matU = randn(sizeX,sizeK);
rvecN = [ 51, 53 ];
[ matV, matX, matS ] = spanspace( matU, rvecN );
numPts = size(matX,2);
for n=1:numPts
	rvecZ(n) = sqrt(min([ ...
	  sum((matX(:,n)).^2), ...
	  sum((matX(:,n)-matU(:,1)).^2), ...
	  sum((matX(:,n)-matU(:,2)).^2) ]));
end
gridZ = reshape( rvecZ, rvecN(2), rvecN(1) );
gridX = reshape( matS(1,:), rvecN(2), rvecN(1) );
gridY = reshape( matS(2,:), rvecN(2), rvecN(1) );
contour( gridX, gridY, sqrt(gridZ), 31 );
hold on;
vecU1 = matU(:,1);
vecU2 = matU(:,2);
vecV1 = matV(:,1);
vecV2 = matV(:,2);
plot( [ 0.0, vecV1'*vecU1 ], [ 0.0, vecV2'*vecU1 ], 'ko-' );
plot( [ 0.0, vecV1'*vecU2 ], [ 0.0, vecV2'*vecU2 ], 'kx-' );
%axis equal;
grid on;
hold off;
toc;
