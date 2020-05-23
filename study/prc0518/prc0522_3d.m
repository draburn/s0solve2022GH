myclear;
%
numXVals = 101;
numYVals = 101;
numZVals = 5;
%
sizeX = 3;
sizeF = 3;
%
if (0)
	funcPrm.matM = eye(sizeF,sizeX);
	funcPrm.vecV = [1;0;0];
	funcPrm.vecA = [2;0;0];
	funcPrm.cQuad = -1.5;
	vecXRoot = [-0.2;0;0];
	ax = [-0.3,1.2,-0.5,0.5];
	z0 = -0.5;
	z1 = 0.5;
	zVals = linspace(z0,z1,numZVals);
	xShiftVals = linspace(0.0,0.0,numZVals);
	yShiftVals = linspace(0.0,0.0,numZVals);
	vecXM = [ 0.789; 0.0; 0.0 ]; % From manual search.
elseif (0)
	funcPrm.vecV = [1;0;0] + 0.5*randn(sizeX,1);
	funcPrm.vecV /= norm(funcPrm.vecV);
	funcPrm.matM = eye(sizeF,sizeX) + 0.3*randn(sizeF,sizeX);
	funcPrm.matM = funcPrm.matM * ( eye(sizeX,sizeX) - 0.5*(funcPrm.vecV*(funcPrm.vecV')) );
	funcPrm.vecA = [2;0;0] + 0.5*randn(sizeF,1);
	funcPrm.cQuad = -2;
	vecXRoot = [-0.2;0;0];
	ax = [-0.5,2.5,-1,1];
	z0 = -0.5; z1 = 0.5;
	%
	% Three roots...
	%ax = [2.005,2.007,0.219,0.221];
	%z0 = -0.3183; z1 = -0.3184;
	vecXM = [ 2.006; 0.220; -0.318 ];
	%
	% Another root!
	%ax = [0.661,0.663,0.085,0.087];
	%z0 = -0.1245; z1 = -0.1243;
	%vecXM = [ 0.662; 0.086; -0.124 ];
	%
	%
	zVals = linspace(z0,z1,numZVals);
	xShiftVals = linspace(0.0,0.0,numZVals);
	yShiftVals = linspace(0.0,0.0,numZVals);
else
	funcPrm.matM = eye(sizeF,sizeX) + 0.2*randn(sizeF,sizeX);
	funcPrm.vecV = [1;0;0] + 0.2*randn(sizeX,1);
	funcPrm.vecA = [2;0;0] + 0.2*randn(sizeF,1);
	funcPrm.cQuad = -1.3;
	vecXRoot = [-0.2;0;0];
	%
	% Good for FN;
	ax = [-0.3,1.2,-0.5,0.5];
	z0 = 0.0; z1 = 0.140;
	%
	% Good for PFN
	ax = [-0.3,0.9,-0.2,0.1];
	z0 = 0.0; z1 = 0.140;
	%
	%ax = [0.672,0.674,-0.147,-0.146];
	%z0 = 0.139; z1 = 0.141;
	vecXM = [ 0.673; -0.146; 0.140 ]; % From manual search.
	%vecXM = [ 0.2; -0.3; 0.140 ]; % An arbitrary value.
	%
	zVals = linspace(z0,z1,numZVals);
	xShiftVals = linspace(0.0,0.0,numZVals);
	yShiftVals = linspace(0.0,0.0,numZVals);
end
vecFPreRoot = oneComponentCubic(vecXRoot,funcPrm);
funchF = @(x)( oneComponentCubic(x,funcPrm) - repmat(vecFPreRoot,[1,size(x,2)]) );
%%%funchF = @(x)( x );
vecFM = funchF(vecXM);
vecFMHat = vecFM/norm(vecFM);
assert( abs(norm(vecFMHat)-1.0) < sqrt(eps) );
matP = eye(sizeF,sizeF) - (vecFMHat*(vecFMHat'));
%
xVals0 = linspace(ax(1),ax(2),numXVals);
yVals0 = linspace(ax(3),ax(4),numYVals);
[ gridX0, gridY0 ] = ndgrid( xVals0, yVals0 );
for n=1:max(size(zVals));
	z = zVals(n);
	gridX = gridX0 + xShiftVals(n);
	gridY = gridY0 + yShiftVals(n);
	rvecX = reshape( gridX, [1,numXVals*numYVals] );
	rvecY = reshape( gridY, [1,numXVals*numYVals] );
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
	dat(n).gridX = gridX;
	dat(n).gridY = gridY;
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
		pfLo = mymin( gridPFN );
		pfHi = mymax( gridPFN );
	else
		fLoSigned = mymin( fLoSigned, gridF1, gridF2, gridF3, gridFN, gridPFN );
		fLo = mymin( fLo, gridFN );
		fHi = mymax( fHi, gridFN );
		pfLo = mymin( pfLo, gridPFN );
		pfHi = mymax( pfHi, gridPFN );
	end
	%
	gridDFNDX = (gridFN(3:end,2:end-1)-gridFN(1:end-2,2:end-1)) ...
	  ./(gridX(3:end,2:end-1)-gridX(1:end-2,2:end-1));
	gridDFNDY = (gridFN(2:end-1,3:end)-gridFN(2:end-1,1:end-2)) ...
	  ./(gridY(2:end-1,3:end)-gridY(2:end-1,1:end-2));
	gridGNC = sqrt( gridDFNDX.^2 + gridDFNDY.^2 );
	gridXC = gridX(2:end-1,2:end-1);
	gridYC = gridY(2:end-1,2:end-1);
	dat(n).gridGNC = gridGNC;
	dat(n).gridXC = gridXC;
	dat(n).gridYC = gridYC;
	if (1==n)
		gLo = mymin(gridGNC);
		gHi = mymax(gridGNC);
	else
		gLo = mymin(gLo,gridGNC);
		gHi = mymax(gHi,gridGNC);
	end
end
%
if (0) % Viz FN on each slice...
for n=1:max(size(zVals));
	figIndex++; figure(figIndex);
	mycontour( dat(n).gridX, dat(n).gridY, dat(n).gridFN, fLo, fHi, mycmap, @(f)asinh(10*f)/10, 30 );
	title(sprintf("z = %f",dat(n).z));
end
return;
end
%
if (0) % Viz GN on each slice...
for n=1:max(size(zVals));
	figIndex++; figure(figIndex);
	funchViz = @(g)asinh(10*g)/10;
	%funchViz = @(g)g;
	mycontour( dat(n).gridXC, dat(n).gridYC, dat(n).gridGNC, gLo, gHi, mycmap, funchViz, 30 );
	title(sprintf("z = %f",dat(n).z));
end
return;
end
%
%
if (1) % Viz PFN on each slice...
for n=1:max(size(zVals));
	figIndex++; figure(figIndex);
	mycontour( dat(n).gridX, dat(n).gridY, dat(n).gridPFN, pfLo, pfHi, mycmap(12), @(f)sqrt(f), 20 );
	title(sprintf("z = %f",dat(n).z));
end
return;
end
