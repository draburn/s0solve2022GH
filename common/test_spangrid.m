clear;
tic();
numFigs = 0;
%
matV = eye(2,2);
matSLoHi = [0,2;1,3];
rvecN = [10,20];
%
aryX = spangrid( matV, matSLoHi, rvecN );
%
numFigs++; figure(numFigs);
imagesc( reshape(aryX(1,:,:),rvecN) );
axis equal;
grid on;
%
numFigs++; figure(numFigs);
imagesc( reshape(aryX(2,:,:),rvecN) );
axis equal;
grid on;
%
%
%
matSLoHi = [0,2,4;1,3,5];
%matV = eye(3,3);
matV = myorth(randn(5,3));
rvecN = [10,20,25];
%
aryX = spangrid( matV, matSLoHi, rvecN );
%
numFigs++; figure(numFigs);
plot( ...
  reshape(aryX(1,:,1,1),1,[]), 'o-', ...
  reshape(aryX(2,1,:,1),1,[]), 'o-', ...
  reshape(aryX(3,1,1,:),1,[]), 'o-' );
grid on;
%
numFigs++; figure(numFigs);
plot( ...
  reshape(aryX(1,:,1,1),1,[]), 'o-', ...
  reshape(aryX(2,:,1,1),1,[]), 'o-', ...
  reshape(aryX(3,:,1,1),1,[]), 'o-', ...
  reshape(aryX(4,:,1,1),1,[]), 'o-', ...
  reshape(aryX(5,:,1,1),1,[]), 'o-' );
grid on;
%
numFigs++; figure(numFigs);
imagesc( reshape(aryX(1,:,:,1),[rvecN(1),rvecN(2)]) );
axis equal;
grid on;
%
numFigs++; figure(numFigs);
imagesc( reshape(aryX(1,:,1,:),[rvecN(1),rvecN(3)]) );
axis equal;
grid on;
%
numFigs++; figure(numFigs);
imagesc( reshape(aryX(1,1,:,:),[rvecN(2),rvecN(3)]) );
axis equal;
grid on;
%
numFigs++; figure(numFigs);
imagesc( reshape(aryX(2,1,:,:),[rvecN(2),rvecN(3)]) );
axis equal;
grid on;
%
numFigs++; figure(numFigs);
imagesc( reshape(aryX(3,1,:,:),[rvecN(2),rvecN(3)]) );
axis equal;
grid on;
%
return;



%
sizeX = 5;
sizeK = 4;
randnState = mod(round(1E6*time),1E6);
%randnState = 0;
%randnState = 470733;
echo_randnState = randnState
randn( "state", randnState );
matU = randn(sizeX,sizeK);
rvecN = [ 15, 17, 19, 7 ];
%
[ matV, matX, matS ] = spanspace( matU, rvecN );
numPts = rvecN(1)*rvecN(2);
for n=1:numPts
	rvecZ(n) = sqrt(min([ ...
	  sum((matX(:,n)).^2), ...
	  sum((matX(:,n)-matU(:,1)).^2), ...
	  sum((matX(:,n)-matU(:,2)).^2) ]));
end
matS1 = matS(1,1:numPts);
matS2 = matS(2,1:numPts);
gridZ = reshape( rvecZ, rvecN(1), rvecN(2) );
gridX = reshape( matS1, rvecN(1), rvecN(2) );
gridY = reshape( matS2, rvecN(1), rvecN(2) );
%
numFigs++; figure(numFigs);
hold off;
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
%
