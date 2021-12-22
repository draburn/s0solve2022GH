myclear;
%setprngstates(64036384);
tic();
%
sizeX = 2
sizeF = 2
numX1Vals = 301;
numX2Vals = 301;
x1Lo = -5.0;
x1Hi =  5.0;
x2Lo = -5.0;
x2Hi =  5.0;
x1Vals = linspace(x1Lo,x1Hi,numX1Vals);
x2Vals = linspace(x2Lo,x2Hi,numX2Vals);
[ gridX1, gridX2 ] = ndgrid( x1Vals, x2Vals );
numVals = numX1Vals*numX2Vals;
matX(1,:) = reshape( gridX1, [1,numVals] );
matX(2,:) = reshape( gridX2, [1,numVals] );
%
matC1 = randn(sizeF,sizeX);
matF = matC1*matX;
%
aryC2 = randn(sizeF,sizeX,sizeX)/2.0;
for n=1:sizeX
for m=1:n
	matF += aryC2(:,n,m) ...
	 * ( matX(n,:) .* matX(m,:) );
end
end
%
aryC3 = randn(sizeF,sizeX,sizeX,sizeX)/6.0;
for n=1:sizeX
for m=1:n
for l=1:m
	matF += aryC3(:,n,m,l) ...
	  * ( matX(n,:) .* matX(m,:) .* matX(l,:) );
end
end
end
%
aryC4 = randn(sizeF,sizeX,sizeX,sizeX,sizeX)/24.0;
for n=1:sizeX
for m=1:n
for l=1:m
for k=1:l
	matF += aryC4(:,n,m,l,k) ...
	  * ( matX(n,:) .* matX(m,:) .* matX(l,:) .* matX(k,:));
end
end
end
end
%
aryC5 = randn(sizeF,sizeX,sizeX,sizeX,sizeX,sizeX)/120.0;
for n=1:sizeX
for m=1:n
for l=1:m
for k=1:l
for j=1:k
	matF += aryC5(:,n,m,l,k,j) ...
	  * ( matX(n,:) .* matX(m,:) .* matX(l,:) .* matX(k,:) .* matX(j,:) );
end
end
end
end
end
%
toc;
tic;
%
gridF1 = reshape( matF(1,:), [numX1Vals,numX2Vals] );
gridF2 = reshape( matF(2,:), [numX1Vals,numX2Vals] );
gridFNorm = sqrt( (gridF1.^2) + (gridF2.^2) );
%
%funchViz = @(f)( f );
%funchViz = @(f)( asinh(f*10.0)/10.0 );
funchViz = @(f)( asinh(f) );
cMap = 0.5 + (0.5*jet(1000));
cMap(end,:) *= 0.4;
cMap(end,:) += 0.6;
cMap(1,:) *= 0.4;
numContours = 20;
%
figIndex++; figure(figIndex);
contourf(gridX1,gridX2,funchViz(gridF1),numContours);
axis equal;
colormap(cMap);
xlabel("x1");
ylabel("x2");
title("funchViz(F1)");
grid on;
%
figIndex++; figure(figIndex);
contourf(gridX1,gridX2,funchViz(gridF2),numContours);
axis equal;
colormap(cMap);
xlabel("x1");
ylabel("x2");
title("funchViz(F2)");
grid on;
%
figIndex++; figure(figIndex);
contourf(gridX1,gridX2,funchViz(gridFNorm),numContours);
axis equal;
colormap(cMap);
xlabel("x1");
ylabel("x2");
title("funchViz(FNorm)");
grid on;
toc;
