clear;
thisFile = "x3_muEscape";
commondefs;
numFigs = 0;
startTime = time();
%
numPts = 101;
doMuVariationPlots = true;
switch (200)
case 0
	x = linspace( 0.0, 4.0, numPts );
	%mu = [ 0.0, 0.5, 3.133, 3.1334, 5.0 ];
	%mu = [ 0.0, 0.5, 3.13, 3.13325, 5.0 ];
	%%%mu = linspace( 3.1332598, 3.1332599, 11 );
	mu = linspace( 3.1332, 3.1333, 11 );
	f = 20 - x - x.^3 + 0.2*x.^4;
case 1
	%x = linspace( -3.0, 6.0, numPts );
	x = linspace( 0.0, 6.0, numPts );
	mu = linspace( 0.0, 10.0, 11 );
	f = 20 - x - x.^3 + 0.2*x.^4;
case 2
	x = linspace( 0.5, 2.5, numPts );
	mu = linspace( 3.133, 3.134, 11 );
	f = 20 - x - x.^3 + 0.2*x.^4;
case 3
	x = linspace( 0.0, 1.5, numPts );
	mu = linspace( 3.133, 3.134, 2 );
	f = 20 - x - x.^3 + 0.2*x.^4;
case 10
	x = linspace( 0.0, 5.0, numPts );
	warning( "This is a different f!" );
	f = 1 - x - x.^3 + 0.25*x.^4;
	%mu = linspace( 0.0, 10.0, 11 );
	mu = linspace( 2.999, 3.0, 11 );
case 100
	x = linspace( -1.0, 1.0, numPts );
	mu = linspace( 0.0, 10.0, 11 );
	f = ( 1 + x - 2*(x.^2) ).^2;
case 200
	% Skips
	x = linspace( 0.0, 10.0, numPts );
	mu = linspace( 0.0, 10.0, 11 );
	x0 = 5.0;
	c3 = 1.0; % impacts left/right shift of f', right?
	c1 = -10.0; % up/down shift of f'.
	alpha = 5.0; % amplitude of f' non-monotonicity.
	x0mu = 0.0;
	switch (1)
	case 0
		% Baseline. Skips.
	case 1
		c1 = -50.0 % No skip.
	case 2
		c1 = -50.0
		x0mu = -10.0
	otherwise
		error( "Invalid sub-case." );
	end
	c2 = ((c3^2)/3.0) - alpha;
	c0 = 200.0; % Doesn't impact f', but could increase to keep f >= 0.
	y = x - x0;
	f = ((y.^4)/4.0) + (c3*(y.^3)/3.0) + (c2*(y.^2)/2.0) + (c1*y) + c0;
otherwise
	error( "Invalid case." );
end
%
%%%g = 0.5*(x.^2);
g = 0.5*((x-x0mu).^2);
numCurves = max(size(mu));
curveColors = 0.6*jet(numCurves);
curveSymbols = [ 'o', 'x', '^', 'v', '*', 's', 'p' ];
for n=1:numCurves
	curveDat(n).x = x;
	h = f+mu(n)*g;
	curveDat(n).h = h;
	curveDat(n).dhdx = diff(h)./diff(x);
	curveDat(n).centx = cent(x);
	%
	iOfFirstMin = 1;
	while (iOfFirstMin<max(size(x)))
		if ( h(iOfFirstMin+1) > h(iOfFirstMin) )
			break;
		end
		iOfFirstMin++;
	end
	curveDat(n).iOfFirstMin = iOfFirstMin;
	curveDat(n).xOfFirstMin = x(iOfFirstMin);
	curveDat(n).hOfFirstMin = h(iOfFirstMin);
	%
	curveDat(n).plot_color = curveColors(n,:);
	curveDat(n).plot_markerSize = 1;
	curveDat(n).plot_markerStyle = curveSymbols(1+mod(n-1,max(size(curveSymbols))));
	curveDat(n).plot_lineWidth = 1;
	curveDat(n).plot_lineStyle = '-';
	curveDat(n).plot_bigMarkerSize = 16;
	curveDat(n).plot_bigLineWidth = 4;
	if ( 0.0 == mu(n) )
		curveDat(n).plot_color = [ 0.0, 0.0, 0.0 ];
	curveDat(n).plot_lineWidth = 3;
	end
	%
	curveDat(n).strName = sprintf("%0.3e",mu(n));
end
%
if (doMuVariationPlots)
numFigs++; figure(numFigs);
n = 1;
plot( ...
  curveDat(n).xOfFirstMin, curveDat(n).hOfFirstMin, ...
  curveDat(n).plot_markerStyle, ...
  'linewidth', curveDat(n).plot_bigLineWidth, ...
  'color', curveDat(n).plot_color, ...
  'markersize', curveDat(n).plot_bigMarkerSize );
hold on;
cellAry_legend(1,n) = { curveDat(n).strName };
for n=2:numCurves
	 plot( ...
	  curveDat(n).xOfFirstMin, curveDat(n).hOfFirstMin, ...
	  curveDat(n).plot_markerStyle, ...
	  'linewidth', curveDat(n).plot_bigLineWidth, ...
	  'color', curveDat(n).plot_color, ...
	  'markersize', curveDat(n).plot_bigMarkerSize );
	cellAry_legend(1,n) = { curveDat(n).strName };
end
grid on;
legend( cellAry_legend, "location", "northwestoutside" );
n = 1;
plot( ...
  curveDat(n).x, curveDat(n).h, ...
  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
  'linewidth', curveDat(n).plot_lineWidth, ...
  'color', curveDat(n).plot_color, ...
  'markersize', curveDat(n).plot_markerSize );
for n=2:numCurves
	plot( ...
	  curveDat(n).x, curveDat(n).h, ...
	  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
	  'linewidth', curveDat(n).plot_lineWidth, ...
	  'color', curveDat(n).plot_color, ...
	  'markersize', curveDat(n).plot_markerSize );
end
hold off;
xlabel( "x" );
ylabel( "omega" );
title( "omega vs x" );
axis([0,10,0,1200]);
%
numFigs++; figure(numFigs);
n = 1;
plot( ...
  curveDat(n).xOfFirstMin, 0.0, ...
  curveDat(n).plot_markerStyle, ...
  'linewidth', curveDat(n).plot_bigLineWidth, ...
  'color', curveDat(n).plot_color, ...
  'markersize', curveDat(n).plot_bigMarkerSize );
hold on;
cellAry_legend(1,n) = { curveDat(n).strName };
for n=2:numCurves
	 plot( ...
	  curveDat(n).xOfFirstMin, 0.0, ...
	  curveDat(n).plot_markerStyle, ...
	  'linewidth', curveDat(n).plot_bigLineWidth, ...
	  'color', curveDat(n).plot_color, ...
	  'markersize', curveDat(n).plot_bigMarkerSize );
	cellAry_legend(1,n) = { curveDat(n).strName };
end
grid on;
legend( cellAry_legend, "location", "northwestoutside" );
n = 1;
plot( ...
  curveDat(n).centx, curveDat(n).dhdx, ...
  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
  'linewidth', curveDat(n).plot_lineWidth, ...
  'color', curveDat(n).plot_color, ...
  'markersize', curveDat(n).plot_markerSize );
for n=2:numCurves
	plot( ...
	  curveDat(n).centx, curveDat(n).dhdx, ...
	  [ curveDat(n).plot_markerStyle, curveDat(n).plot_lineStyle], ...
	  'linewidth', curveDat(n).plot_lineWidth, ...
	  'color', curveDat(n).plot_color, ...
	  'markersize', curveDat(n).plot_markerSize );
end
hold off;
xlabel( "x" );
ylabel( "d/dx omega" );
title( "d/dx omega vs x" );


if (0)
numFigs++; figure(numFigs);
plot( cent(x), (-diff(f)./diff(x))./cent(x), 'o-' );
grid on;
xlabel( "x" );
ylabel( "-g / x" );
title( "-g / x vs x" );
end

end


if (0)
numFigs++; figure(numFigs);
plot( ...
  cent(x), -diff(f)./diff(x), 'o-', ...
  x, 3.14*x, '-', 'linewidth', 2 );
grid on;
xlabel( "x" );
ylabel( "-g" );
title( "-g vs x" );
end


df = diff(f)./diff(x);
cx = cent(x);
ddf = diff(df)./diff(cx);
ccx = cent(cent(x));
%
numFigs++; figure(numFigs);
plot( ...
  x, f, 'o-' );
grid on;
xlabel( "x" );
ylabel( "omega" );
title( "omega vs x" );
%
if (0)
df = diff(f)./diff(x);
cx = cent(x);
ddf = diff(df)./diff(cx);
ccx = cent(cent(x));
numFigs++; figure(numFigs);
plot( ...
  ccx, ccx.*(-ddf), 'o-', ...
  cx, -df, 'x-' );
grid on;
xlabel( "x" );
ylabel( "x*(d/dx(-g)) and -g" );
title( "-g vs x" );
else
	numFigs++; figure(numFigs);
	plot( cx, df, 'o-' );
	grid on;
	xlabel( "x" );
	ylabel( "omega'" );
	title( "omega' vs x" );
end

if (0)
numFigs++; figure(numFigs);
plot( cent(x), (-diff(f)./diff(x))./(cent(x).^3), 'o-' );
grid on;
xlabel( "x" );
ylabel( "-(d/dx omega ) / x\^3" );
title( "-(d/dx omega ) / x\^3 vs x" );
end
