myclear;
%
numPts = 0;
pts = [];
%numPts++; pts(numPts,:) = [ -0.7, 1.3, 0.4, 0.0, -2.0 ];
%numPts++; pts(numPts,:) = [  0.7, 1.3, 0.4, 0.0,  8.0 ];
%numPts++; pts(numPts,:) = [ -0.8, 1.5, 0.4, -1.0, -9.0 ];
%numPts++; pts(numPts,:) = [  0.8, 1.9, 0.4, -1.0, -7.0 ];
%numPts++; pts(numPts,:) = [ 0.0, 2.0, 0.4, 0.0, -5.0 ];
%numPts++; pts(numPts,:) = [ -2.0, 0.0, 0.3, -10.0, 0.0 ];
%numPts++; pts(numPts,:) = [  0.0, 0.0, 0.3, 0.0, -5.0 ];
%numPts++; pts(numPts,:) = [  2.0, 0.0, 1.0, 2.0, 0.0 ];
%numPts++; pts(numPts,:) = [  1.0, 0.0, 0.5, 20.0, -35 ];
%
numPts++; pts(numPts,:) = [  0.5, 0.0, 0.5, 20.0, 0.0 ];
numPts++; pts(numPts,:) = [  -2.0, 0.0, 0.5, 20.0, 0.0 ];

%
%funchF1 = @(x,y) (foo0514_vizcrit__f1(x,y,pts) - foo0514_vizcrit__f1(0,0,pts));
%funchF2 = @(x,y) (foo0514_vizcrit__f2(x,y,pts) - foo0514_vizcrit__f2(0,0,pts));
funchF1 = @(x,y) foo0514_vizcrit__f1(x,y,pts);
funchF2 = @(x,y) foo0514_vizcrit__f2(x,y,pts);
%ax = [ -2.5, 2.5, -1, 4.0 ];
%ax = [ -2.0, 2.0, 0.0, 4.0 ];
%ax = [ 1, 1.2, -0.1, 0.1 ];
ax = 5*[-1,1,-1,1];
numXVals = 51;
numYVals = 51;
%
xVals = linspace(ax(1),ax(2),numXVals);
yVals = linspace(ax(3),ax(4),numYVals);
[ gridX, gridY ] = ndgrid( xVals, yVals );
%
gridF1 = funchF1(gridX,gridY);
gridF2 = funchF2(gridX,gridY);
gridFNorm = sqrt(gridF1.^2+gridF2.^2);
%
d_xVals = cent(xVals);
d_yVals = cent(yVals);
d_gridX = gridX(2:end-1,2:end-1);
d_gridY = gridY(2:end-1,2:end-1);
dx_gridFNorm = gridFNorm(3:end,2:end-1)-gridFNorm(1:end-2,2:end-1);
dy_gridFNorm = gridFNorm(2:end-1,3:end)-gridFNorm(2:end-1,1:end-2);
d_gridFNorm = sqrt( dx_gridFNorm.^2 + dy_gridFNorm.^2 );
%
%dx_gridFNorm= gridFNorm(2:end,:)-gridFNorm(1:end-1,:);
%d_gridFNorm_dxFlag = double(dx_gridFNorm(2:end,2:end-1).*dx_gridFNorm(1:end-1,2:end-1) <= 0.0);
%dy_gridFNorm= gridFNorm(:,2:end)-gridFNorm(:,1:end-1);
%d_gridFNorm_dyFlag = double(dy_gridFNorm(2:end-1,2:end).*dy_gridFNorm(2:end-1,1:end-1) <= 0.0);
%
foo0514_vizcrit__viz;
