clear;
%setprngstates
setprngstates(0);
numFigs = 0;
thisFile = "chaoticSignal0712"
%
r = @(n) randn(n,1);
s = @(x)( x + sin(x) );
t = @(x,c)( s( c(1) + abs(c(2))*s( c(3) + abs(c(4))*s( c(5) + abs(c(6))*x ) ) ) );
u = @(x,c)( t( x, c(1:6) ) + t( x, c(7:12) ) );
v = @(x,c)( u(u( x, c(1:12) ),c(13:24)) );
y = @(x,c)( v(x,c) );
%xzLo = -1e8;
%xzHi = 1e8;
%z = @(x,c)( y(x,c) * (xzHi-xzLo)/( y(xzHi,c) - y(xzLo,c) ) );
%
x = 10002+linspace(-1,1,10000)*3;
%plot( x, y(x,2*pi*r(24)), 'o-' );
plot( x, cos(2*pi*r(1)+y(x,2*pi*r(24))), "o-" );
%plot( x, cos(p(1)+z(x,4.0*p(24))), 'o-' );
grid on;
return;
