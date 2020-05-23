myclear;
%
numXVals = 51;
numYVals = 51;
numZVals = 5;
%
sizeX = 3;
sizeF = 3;
%
if (1)
	funcPrm.matM = eye(sizeF,sizeX);
	funcPrm.vecV = [1;0;0];
	funcPrm.vecA = [2;0;0];
	funcPrm.cQuad = -1.5;
	vecXRoot = [-0.2;0;0];
	ax = [-0.3,1.2,-0.5,0.5];
	z0 = -0.5;
	z1 = 0.5;
	zVals = linspace(z0,z1,numZVals);
	vecXM = [ 0.789; 0.0; 0.0 ]; % From manual search.
end
vecFPreRoot = oneComponentCubic(vecXRoot,funcPrm);
funchF = @(x)( oneComponentCubic(x,funcPrm) - repmat(vecFPreRoot,[1,size(x,2)]) );
%%%funchF = @(x)( x );
vecFM = funchF(vecXM);
vecFMHat = vecFM/norm(vecFM);
assert( abs(norm(vecFMHat)-1.0) < sqrt(eps) );
matP = eye(sizeF,sizeF) - (vecFMHat*(vecFMHat'));
%
xVals = linspace(ax(1),ax(2),numXVals);
yVals = linspace(ax(3),ax(4),numYVals);
[ gridX, gridY ] = ndgrid( xVals, yVals );
rvecX = reshape( gridX, [1,numXVals*numYVals] );
rvecY = reshape( gridY, [1,numXVals*numYVals] );
for n=1:max(size(zVals));
	z = zVals(n);
	%
	matF = funchF([rvecX;rvecY;z+(0*rvecX)]);
	matPF = matP * matF;
	rvecFN = sqrt(sum(matF.^2,1));
	rvecPFN = sqrt(sum(matPF.^2,1));
	gridF1 = reshape( matF(1,:), [numXVals,numYVals] );
	gridF2 = reshape( matF(2,:), [numXVals,numYVals] );
	gridF3 = reshape( matF(3,:), [numXVals,numYVals] );
	gridFN = reshape( rvecFN, [numXVals,numYVals] );
	gridPFN = reshape( rvecPFN, [numXVals,numYVals] );
	%
	dat(n).z = z;
	dat(n).gridF1 = gridF1;
	dat(n).gridF2 = gridF2;
	dat(n).gridF3 = gridF3;
	dat(n).gridFN = gridFN;
	dat(n).gridPFN = gridPFN;
	%
	if (1==n)
		fLoSigned = mymin( gridF1, gridF2, gridF3, gridFN, gridPFN );
		fLo = mymin( gridFN );
		fHi = mymax( gridFN );
	else
		fLoSigned = mymin( fLoSigned, gridF1, gridF2, gridF3, gridFN, gridPFN );
		fLo = mymin( fLo, gridFN );
		fHi = mymax( fHi, gridFN );
	end
end
%
if (1) % Viz FN on each slice...
for n=1:max(size(zVals));
	figIndex++; figure(figIndex);
	mycontour( gridX, gridY, dat(n).gridFN, fLo, fHi, mycmap, @(z)(z), 30 );
end
return;
end
%
%
if (1) % Viz PFN on each slice...
for n=1:max(size(zVals));
	figIndex++; figure(figIndex);
	mycontour( gridX, gridY, dat(n).gridPFN, fLo, fHi, mycmap, @(z)(z), 30 );
end
return;
end
