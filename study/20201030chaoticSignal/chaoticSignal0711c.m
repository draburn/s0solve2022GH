clear;
setprngstates
%setprngstates(0);
numFigs = 0;
thisFile = "chaoticSignal0711c"
%
r = @(n) randn(n,1);
s = @(x)( x + sin(x) );
%t = @(x,c)( s( c(1) + c(2)*s( c(3) + c(4)*s( c(5) + c(6)*x ) ) ) );
t = @(x,c)( s( c(1) + abs(c(2))*s( c(3) + abs(c(4))*s( c(5) + abs(c(6))*x ) ) ) );
u = @(x,c)( t( x, c(1:6) ) + t( x, c(7:12) ) );
v = @(x,c)( u(u( x, c(1:12) ),c(13:24)) );
f = @(x,c)( cos(c(1)+v(x,c(2:25))) );
a = @(x,c)( 1.0 + f(x,c) );
%
x = 10002+linspace(-1,1,10000)*1;
%plot( x, a(x/10.0,2*pi*r(25)).*f(x,2*pi*r(25)), "o-" );
plot( x, f(x,2*pi*r(25)), "o-" );
grid on;
return;
