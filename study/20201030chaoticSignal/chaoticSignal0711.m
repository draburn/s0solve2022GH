clear;
setprngstates
setprngstates(0);
numFigs = 0;


%d = @(x)( x - round(x) );
%p = @(x)( sign(d(x)) .* abs(2.0*d(x)).^4 );
%y = @(x)( p(x) + 2.0*round(x) );
r = @(n) randn(n,1);
s = @(x)( x + sin(x) );
%t = @(x,c)( s( c(1) + c(2)*s( c(3) + c(4)*s( c(5) + c(6)*x ) ) ) );
t = @(x,c)( s( c(1) + abs(c(2))*s( c(3) + abs(c(4))*s( c(5) + abs(c(6))*x ) ) ) );
u = @(x,c)( t( x, c(1:6) ) + t( x, c(7:12) ) );
v = @(x,c)( u(u( x, c(1:12) ),c(13:24)) );
y = @(x,c)( v(x,c) );
%xzLo = -1e8;
%xzHi = 1e8;
%z = @(x,c)( y(x,c) * (xzHi-xzLo)/( y(xzHi,c) - y(xzLo,c) ) );
%
x = linspace(-1,1,10000)*10;
%plot( x, y(x,10*r(24)), 'o-' );
plot( x, cos(y(x,4*r(24))), 'o-' );
grid on;
return;



%z = @(x)( s(s(x)) );
z = @(x)( s(s(s(x))) );
%z = @(x)( s(s(s(s(s(x))))) );
%plot( x, cos(x+10*z(x/3)), 'o-' );
plot( x, cos(2+0.2*x+10*z(0.2+x/3)).*cos(3+0.3*x+4*z(0.6+x/4.0)), 'o-' );
%plot( x, (1.0+sin(4+x+sin(x+5))).^2 .* cos(2+0.2*x+10*z(3+x/3)), 'o-' );

%plot(x, cos(3*(x+z(x))), 'o-' );

%%%%%%%%
%plot( x, cos(x+z(10*x)), 'o-' );
%plot( x, s(s(s(x))), 'o-' );
%plot( x, cos(2*z(x)+z(1.2*x)+z(1.3*x)), 'o-' );
%plot( x, s(0.1+0.5*s(0.2+x/3.0)+0.6*s(0.3+x/4.0)), 'o-' );
%plot( x, y(0.3+3*y(x/3)), 'o-' );
%plot(x, y(x*10.0/11.0)+y(x*10.0/13.0), 'o-' );
%plot( x, cos(1+50*y(a(x))), 'o-' );
grid on;
return;


r = @(n) randn(n,1);
y = @(x)( x + sin(x) );
x = linspace(-10,10,1000);
plot( x, cos(1+10*y(x)+y(x/10)), 'o-' );
grid on;
return;
y1 = @(x1,x2,c4)( 0 );
y2 = @(x1,x2,x3,x4,c12)(  y1( ...
  y1(x1,x2,c12(1:4)), ...
  y1(x3,x4,c12(5:8)), ...
  c12(9:12) )  );
y1x = @(x,c4) y1(x,x,c4);
y2x = @(x,c12) y2(x,x,x,x,c12);
plot( x, y2x(x,r(12)), 'o-' );
grid on;
return;

omega = (1.0 + cos(x/13.0)).^2;
a = (1.0+cos(x/11.0)).^2;
f = a.*cos( omega.*x );
numFigs++; figure(numFigs);
plot( x, f, 'o-' );
grid on;

return;
%

fa1 = @(x)( randn*x + randn*cos(randn+randn*x) );
fa2 = @(x)( randn*x + randn*cos(randn+randn*x) );
fa3 = @(x)( randn*x + randn*cos(randn+randn*x) );
fa4 = @(x)( randn*x + randn*cos(randn+randn*x) );
fa5 = @(x)( randn*x + randn*cos(randn+randn*x) );
fa6 = @(x)( randn*x + randn*cos(randn+randn*x) );
fa7 = @(x)( randn*x + randn*cos(randn+randn*x) );
fa8 = @(x)( randn*x + randn*cos(randn+randn*x) );
fb1 = @(x)( randn*fa1(x) + randn*cos(randn+randn*fa2(x)) );
fb2 = @(x)( randn*fa3(x) + randn*cos(randn+randn*fa4(x)) );
fb3 = @(x)( randn*fa5(x) + randn*cos(randn+randn*fa6(x)) );
fb4 = @(x)( randn*fa7(x) + randn*cos(randn+randn*fa8(x)) );
%
%fc1 = @(x)( (2.0+cos(fb1(x))) .* cos( (1.0+cos(fb2(x/10.0))) .* x ) );

%
x = linspace(-100,100,1000);
numFigs++; figure(numFigs);
plot( x, f, 'o-' );
grid on;

return;
numYTerms = 4;
yScale = 50.0;
ya0_orig = 0.0;
ya_orig = abs(randn(1,numYTerms));
ya0 = ya0_orig/(ya0_orig+sum(ya_orig));
ya = ya_orig/(ya0_orig+sum(ya_orig));
yf = 2*pi*randn(1,numYTerms)/yScale;
yd = 2*pi*rand(1,numYTerms);
funch_y = @(x)( ya0*x + ...
 ya(1)*( x + cos(yf(1)*x+yd(1))/yf(1) ) + ...
 ya(2)*( x + cos(yf(2)*x+yd(2))/yf(2) ) + ...
 ya(3)*( x + cos(yf(3)*x+yd(3))/yf(3) ) + ...
 ya(4)*( x + cos(yf(4)*x+yd(4))/yf(4) ) );
%
numFTerms = 5;
fScale = 1.0;
fa = randn(1,numFTerms);
ff = 2*pi*randn(1,numFTerms)/fScale;
fd = 2*pi*rand(1,numFTerms);
funch_f = @(x)( ...
  fa(1)*cos(ff(1)*x+fd(1)) + ...
  fa(2)*cos(ff(2)*x+fd(2)) + ...
  fa(3)*cos(ff(3)*x+fd(3)) + ...
  fa(4)*cos(ff(4)*x+fd(4)) + ...
  fa(5)*cos(ff(5)*x+fd(5)) );
%
numFigs = 0;
%
xVals = linspace(-500,500,10001);
numFigs++; figure(numFigs);
plot( cent(xVals), diff(funch_y(xVals))./diff(xVals), 'o-' );
grid on;
numFigs++; figure(numFigs);
plot( xVals, xVals, 'kx-', xVals, funch_y(xVals), 'mx-' );
grid on;
numFigs++; figure(numFigs);
plot( xVals, funch_f(xVals), 'gs-' );
grid on;
numFigs++; figure(numFigs);
plot( xVals, funch_f(funch_y(xVals)), 'r^-' );
grid on;
%
yVals = funch_y(xVals);
dydxVals = diff(yVals)./diff(xVals);
[ foo, iOfMin ] = min(dydxVals);
xOfMin = (xVals(iOfMin)+xVals(iOfMin+1))/2.0;
[ foo, iOfMax ] = max(dydxVals);
xOfMax = (xVals(iOfMax)+xVals(iOfMax+1))/2.0;
%
xLo = xOfMin-0.05*(xOfMax-xOfMin);
xHi = xOfMin+0.05*(xOfMax-xOfMin);
xVals = linspace(xLo,xHi,10001);
numFigs++; figure(numFigs);
plot( xVals, funch_f(xVals), 'gs-' );
grid on;
numFigs++; figure(numFigs);
plot( xVals, funch_f(funch_y(xVals)), 'r^-' );
grid on;
%
xLo = xOfMax-0.05*(xOfMax-xOfMin);
xHi = xOfMax+0.05*(xOfMax-xOfMin);
xVals = linspace(xLo,xHi,10001);
numFigs++; figure(numFigs);
plot( xVals, funch_f(xVals), 'gs-' );
grid on;
numFigs++; figure(numFigs);
plot( xVals, funch_f(funch_y(xVals)), 'r^-' );
grid on;
