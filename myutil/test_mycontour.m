myclear;
thisFile = "test_mycontour";
%
funchZ = @(x,y)( x.^2 + y.^2 - exp( -20.0*((x-1.0).^2+y.^2)) );
[ gridX1, gridY1, gridZ1 ] = gridfunch( funchZ, 2, [-5.0,5.0,-5.0,5.0] );
[ gridX2, gridY2, gridZ2 ] = gridfunch( funchZ, 2, [-1.5,1.5,-1,1], 101, 101 );
funchViz = @(z)(asinh(z*10.0)/10.0);
%
figIndex++; figure(figIndex);
contourf( gridX1, gridY1, funchViz(gridZ1), 20 );
axis equal;
colormap(mycmap);
grid on;
%
figIndex++; figure(figIndex);
contourf( gridX2, gridY2, funchViz(gridZ2), 20 );
%axis equal;
colormap(mycmap);
grid on;
%
%
zLo = mymin(gridZ1,gridZ2)
zHi = mymax(gridZ1,gridZ2)
%
figIndex++; figure(figIndex);
mycontour( gridX1, gridY1, gridZ1, zLo, zHi, mycmap, funchViz, 20 );
axis equal;
figIndex++; figure(figIndex);
mycontour( gridX2, gridY2, gridZ2, zLo, zHi, mycmap, funchViz, 20 );
