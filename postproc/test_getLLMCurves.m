clear;
commondefs;
getLLMCurves_setCnsts;
thisFile = "test_getLLMCurvs";
%
vecF = [-1;-1];
matJ = [-1,1;1,10];
matV = eye(2,2);
matW = matJ*matV;
numPts = 50;
curveTypes = [ ...
  GETCURVES_CURVETYPE__NEWTON, ...
  GETCURVES_CURVETYPE__PICARD, ...
  GETCURVES_CURVETYPE__PICARD_SCALED, ...
  GETCURVES_CURVETYPE__GRADDIR, ...
  GETCURVES_CURVETYPE__GRADDIR_SCALED, ...
  GETCURVES_CURVETYPE__LEVCURVE, ...
  GETCURVES_CURVETYPE__LEVCURVE_SCALED ];
%
prm.curveTypes = curveTypes;
[ curveDat, retCode, datOut ] = getLLMCurves( vecF, matV, matW, numPts, prm );
%
numCurves = max(size(curveTypes));
for n=1:numCurves
	matY = curveDat(n).matY;
	numPts = size(matY,2);
	matRes = repmat(vecF,[1,numPts]) + (matW*matY);
	myDat(n).deltaNorm = sqrt(sum(matY.^2,1));
	myDat(n).matRes = matRes;
	myDat(n).omega = 0.5*sum(matRes.^2,1);
end
%
numFigs = 0;
colMap = 0.8*jet(numCurves);
strLegend = [ curveDat(1).strType ];
for n=2:numCurves
	strLegend = [ strLegend; curveDat(n).strType ];
end
%
numFigs++; figure(numFigs);
plot( curveDat(1).matY(1,:), curveDat(1).matY(2,:), ...
  'o-', 'color', colMap(1,:), 'markersize', 3 );
hold on;
for n=2:numCurves
	plot( curveDat(n).matY(1,:), curveDat(n).matY(2,:), ...
	  'o-', 'color', colMap(n,:), 'markersize', 3*n );
end
hold off;
grid on;
legend( strLegend );
%
numFigs++; figure(numFigs);
plot( myDat(1).deltaNorm, myDat(1).omega, ...
  'o-', 'color', colMap(1,:), 'markersize', 3 );
hold on;
for n=2:numCurves
	plot( myDat(n).deltaNorm, myDat(n).omega, ...
	  'o-', 'color', colMap(n,:), 'markersize', 3*n );
end
hold off;
grid on;
legend( strLegend );
