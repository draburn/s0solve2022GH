clear;
tic();
numFigs = 0;
%randnState = mod(round(1E6*time),1E6);
randnState = 0;
%randnState = 470733;
echo_randnState = randnState
randn( "state", randnState );
rand( "state", randnState );
%
%
sizeX = 3;
sizeK = 3;
matU = randn(sizeX,sizeK);
%matU = eye(sizeX,sizeK);
%matU = diag([3,2,1]);
matV = myorth(matU);
matVTU = matV' * matU;
for indexK=1:sizeK
	matSLoHi(1,indexK) = min([ 0.0, matVTU(indexK,:) ]);
	matSLoHi(2,indexK) = max([ 0.0, matVTU(indexK,:) ]);
end
for indexK=1:2
	% Add border for first two dims only.
	sLo = matSLoHi(1,indexK);
	sHi = matSLoHi(2,indexK);
	matSLoHi(1,indexK) = sLo - 0.1*(sHi-sLo);
	matSLoHi(2,indexK) = sHi + 0.1*(sHi-sLo);
end
%rvecN = [ 15, 7, 3 ];
rvecN = [ 51, 51, 5 ];
[ aryX, aryS ] = spangrid( matV, matSLoHi, rvecN );
sizeC = 256;
colmapFull = 0.2+(0.8*jet(sizeC));
%
numPts = prod(rvecN)
matX = reshape(aryX,[sizeX,numPts]);
sizeF = sizeX;
matJ = randn(sizeF,sizeX);
vecX0 = matU*rand(sizeK,1)*0.5;
matF = matJ*( matX - repmat(vecX0,[1,numPts]) );
rvecOmega = sqrt( 0.5*sum(matF.^2,1) );
aryOmega = reshape( rvecOmega, rvecN );
maxOmega = max(max(max(aryOmega)));
minOmega = min(min(min(aryOmega)));
%
for indexZ = 1 : rvecN(3);
	gridX = reshape(aryS(1,:,:,indexZ),rvecN(1:2))';
	gridY = reshape(aryS(2,:,:,indexZ),rvecN(1:2))';
	gridT = reshape(aryOmega(:,:,indexZ),rvecN(1:2))';
	numFigs++; figure(numFigs);
	hold off;
	contourf( gridX, gridY, gridT, 32 );
	minT = min(min(gridT));
	maxT = max(max(gridT));
	cLo = floor(1+((sizeC-1)*(minT-minOmega)/(maxOmega-minOmega)));
	cHi = ceil(1+((sizeC-1)*(maxT-minOmega)/(maxOmega-minOmega)));
	colormap(colmapFull(cLo:cHi,:));
	title(sprintf( "Omega at Z = %f.", aryS(3,1,1,indexZ) ));
	axis equal;
		hold on;
	if (1==indexZ)
		plot( [0.0, matVTU(1,1)], [0.0, matVTU(2,1)], 'ko-', 'linewidth', 3, 'markersize', 10 );
		plot( [0.0, matVTU(1,2)], [0.0, matVTU(2,2)], 'kx-', 'linewidth', 3, 'markersize', 10 );
	else
		plot( [0.0, matVTU(1,1)], [0.0, matVTU(2,1)], 'ko:', 'linewidth', 1, 'markersize', 10 );
		plot( [0.0, matVTU(1,2)], [0.0, matVTU(2,2)], 'kx:', 'linewidth', 1, 'markersize', 10 );
	end
	plot( [0.0, matVTU(1,3)], [0.0, matVTU(2,3)], 'ks:', 'linewidth', 1, 'markersize', 10 );
	sScaled = aryS(3,1,1,indexZ)/aryS(3,1,1,end);
	plot( matVTU(1,3)*sScaled, matVTU(2,3)*sScaled, 'ks', 'linewidth', 3, 'markersize', 12 );
	soluVec = matV'*vecX0;
	plot( soluVec(1), soluVec(2), 'kp', 'linewidth', 1, 'markersize', 15 );
	hold off;
	grid on;
end
%
%
numFigs++; figure(numFigs);
gridX = repmat(linspace(0,1,2),[sizeC,1]);
gridY = repmat(linspace(minOmega,maxOmega,sizeC), [2, 1])';
gridT = repmat(linspace(minOmega,maxOmega,sizeC), [2, 1])';
contourf( gridX, gridY, gridT, 32 );
set(get(gcf,"children"),"xtick",[])
colormap(colmapFull);
title(sprintf( "Color scale (%g to %g)", minOmega, maxOmega ));
grid on;
%
%
return;
%
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
