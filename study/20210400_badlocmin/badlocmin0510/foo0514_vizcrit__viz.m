figIndex = 0;
%
cMap = 0.5+0.5*jet;
cMap(end,:) *= 0.4;
cMap(end,:) += 0.6;
cMap(1,:) *= 0.4;
%
figIndex++; figure(figIndex);
contourf( gridX, gridY, (gridFNorm*100)/100, 20 );
axis equal;
colormap(cMap);
grid on;
return


if (1)
figIndex++; figure(figIndex);
contourf( gridX, gridY, sign(gridF1), 1 );
axis equal;
colormap(0.5+0.5*jet);
grid on;
%
figIndex++; figure(figIndex);
contourf( gridX, gridY, sign(gridF2), 1 );
axis equal;
colormap(0.5+0.5*jet);
grid on;
end
%
if (1)
figIndex++; figure(figIndex);
contourf( d_gridX, d_gridY, sign(dx_gridFNorm), 1 );
axis equal;
colormap(0.5+0.5*jet);
grid on;
%
figIndex++; figure(figIndex);
contourf( d_gridX, d_gridY, sign(dy_gridFNorm), 1 );
axis equal;
colormap(0.5+0.5*jet);
grid on;
end
return
%
figIndex++; figure(figIndex);
%contour( d_gridX, d_gridY, asinh(d_gridFNorm*10)/10, 20 );
contour( d_gridX, d_gridY, d_gridFNorm.^0.3, 20 );
%imagesc( d_xVals, d_yVals, asinh(d_gridFNorm*10)/10 );
axis equal;
colormap(hot(1000));
grid on;
