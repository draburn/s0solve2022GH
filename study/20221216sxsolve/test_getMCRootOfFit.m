clear;
setprngstates();
%
xL = randn()*exp(3.0*randn())
xR = xL + abs(randn()*exp(3.0*randn()))%*(1.0-2.0*(rand()<0.1))
fL = randn()*exp(3.0*randn())
fR = -sign(fL)*abs(randn()*exp(3.0*randn()))%*(1.0-2.0*(rand()<0.1))
dfdxL = randn()*exp(3.0*randn())
dfdxR = sign(dfdxL)*abs(randn()*exp(3.0*randn()))%*(1.0-2.0*(rand()<0.1))
xMCRoot = getMCRootOfFit( xL, fL, dfdxL, xR, fR, dfdxR )
%
sViz = linspace(0.0,1.0,101);
xSpan = xR - xL;
xViz = xL + xSpan*sViz;
foo = [ 3.0, -1.0; -2.0, 1.0 ] * [ fR - fL - dfdxL * xSpan; (dfdxR - dfdxL) * xSpan ];
c0 = fL
c1 = dfdxL * xSpan;
c2 = foo(1);
c3 = foo(2);
fViz = c0 + c1*sViz + c2*(sViz.^2) + c3*(sViz.^3);
fLViz = fL + dfdxL*xSpan*sViz;
fRViz = fR + dfdxR*xSpan*(sViz-1.0);
%
plot( ...
  xViz, fViz, 'o-', ...
  xViz, 0.0*xViz, 'k-', ...
  xL, fL, 'v', 'markerSize', 20, ...
  xR, fR, '^', 'markerSize', 20, ...
  xViz, fLViz, '-', ...
  xViz, fRViz, '-' );
grid on;
if ( isempty(xMCRoot) )
	msg( __FILE__, __LINE__, "MCRoot does not exist." );
else
	msg( __FILE__, __LINE__, sprintf( "MCRoot is at %f.", xMCRoot ) );
	hold on;
	plot( xMCRoot, 0.0, '*', 'markerSize', 20 );
	hold off;
endif

