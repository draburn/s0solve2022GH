myclear;
%
funchA = @(x,y)( (1./(1+(x-1).^2+y.^2)).^8 );
funchB = @(x,y)( (1./(1+(x+1).^2+y.^2)).^8 );
funchC = @(x,y)( (1./(1+(y-1).^2+x.^2)).^8 );
funchD = @(x,y)( (1./(1+(y+1).^2+x.^2)).^8 );
funchF1 = @(x,y)( funchA(x,y) + funchB(x,y) - funchC(x,y) - funchD(x,y));
funchF2 = @(x,y)( 0.01*(x.^2 + y.^2) );
funchOmega = @(x,y)sqrt( funchF1(x,y).^2 + funchF2(x,y).^2 );
%
%funchZ = @(x,y)( funchOmega(x,y) );
vizC = 1E1;
funchZ = @(x,y)( abs(asinh(funchOmega(x,y)*vizC)/vizC) );
figIndex++; figure(figIndex);
%contourfunch(funchZ,2*[-1,1,-1,1],201,201,50)
contourfunch(funchZ,[0,1,0,1],201,201,50)
%contourfunch(funchZ,[-0.1,0.1,-2,2],201,201,50)
%axis equal;
colormap(0.5+0.5*jet(1000));
grid on;
%
%Look at deriv...
