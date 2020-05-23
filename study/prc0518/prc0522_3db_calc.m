numXVals = 201;
numYVals = 201;
numZVals = 5;
minVar = 0.5;
axCoeff = 2.0;
%
if (1)
	% Shows main space btw vecXM and vecXRoot.
	%ax = [-1,1.5,-1,1];
	%z0 = min([vecXRoot(3),vecXM(3)]);
	%z1 = max([vecXRoot(3),vecXM(3)]);
	%
	% Shows PRC from vecXM.
	ax = [-1,0.5,-0.5,1.5];
	z0 = -0.688;
	z1 = max([vecXRoot(3),vecXM(3)]);
	%
	xShift0 = 0.0; xShift1 = 0.0;
	yShift0 = 0.0; yShift1 = 0.0;
else
	temp_xLo = min([vecXRoot(1),vecXM(1)]);
	temp_xHi = max([vecXRoot(1),vecXM(1)]);
	temp_yLo = min([vecXRoot(2),vecXM(2)]);
	temp_yHi = max([vecXRoot(2),vecXM(2)]);
	temp_xMid = (temp_xHi+temp_xLo)/2.0;
	temp_xVar = max([ (temp_xHi-temp_xLo)/2.0, minVar ]);
	temp_yMid = (temp_yHi+temp_yLo)/2.0;
	temp_yVar = max([ (temp_yHi-temp_yLo)/2.0, minVar ]);
	temp_x0 = temp_xMid - axCoeff*temp_xVar;
	temp_x1 = temp_xMid + axCoeff*temp_xVar;
	temp_y0 = temp_yMid - axCoeff*temp_yVar;
	temp_y1 = temp_yMid + axCoeff*temp_yVar;
	ax = [temp_x0,temp_x1,temp_y0,temp_y1]
	%
	z0 = min([vecXRoot(3),vecXM(3)]);
	z1 = max([vecXRoot(3),vecXM(3)]);
	%
	xShift0 = 0.0; xShift1 = 0.0;
	yShift0 = 0.0; yShift1 = 0.0;
end
%
vecFMHat = vecFM/norm(vecFM);
assert( abs(norm(vecFMHat)-1.0) < sqrt(eps) );
matP = eye(sizeF,sizeF) - (vecFMHat*(vecFMHat'));
%
xVals0 = linspace(ax(1),ax(2),numXVals);
yVals0 = linspace(ax(3),ax(4),numYVals);
zVals = linspace(z0,z1,numZVals);
xShiftVals = linspace(xShift0,xShift1,numZVals);
yShiftVals = linspace(yShift0,yShift1,numZVals);
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
