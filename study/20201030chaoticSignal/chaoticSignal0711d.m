clear;
setprngstates
%setprngstates(0);
%setprngstates(78886496);
%setprngstates(20859712);
numFigs = 0;
thisFile = "chaoticSignal0711d";
%
xNormyLo = -100.0;
xNormyHi = 100.0;
rf = @(n) rand(n,1);
rg = @(n) randn(n,1);
s = @(x)( x + sin(x) );
s1u = @(x,rf,rg)( s( rf + rg*s(x) ) );
s1n = @(x,rf,rg)( s1u(x,rf,rg) * (xNormyHi-xNormyLo)/(s1u(xNormyHi,rf,rg)-s1u(xNormyLo,rf,rg) ) );
s1 = @(x,rf,rg) s1n(x,rf,rg);
s2u = @(x,rf,rg)( s1( rf(1) + rg(1)*s1(x,rf(2),rg(2)), rf(3), rg(3) ) );
s2n = @(x,rf,rg)( s2u(x,rf,rg) * (xNormyHi-xNormyLo)/(s2u(xNormyHi,rf,rg)-s2u(xNormyLo,rf,rg) ) );
s2 = @(x,rf,rg) s2n(x,rf,rg);
s3u = @(x,rf,rg)( s2( rf(1) + rg(1)*s2(x,rf(2:4),rg(2:4)), rf(5:7), rg(5:7) ) );
s3n = @(x,rf,rg)( s3u(x,rf,rg) * (xNormyHi-xNormyLo)/(s3u(xNormyHi,rf,rg)-s3u(xNormyLo,rf,rg) ) );
s3 = @(x,rf,rg) s3n(x,rf,rg);

%t = @(x,rf,rg)( s( 2*pi(1) + abs(c(2))*s( c(3) + abs(c(4))*s( c(5) + abs(c(6))*x ) ) ) );
%u = @(x,c)( t( x, c(1:6) ) + t( x, c(7:12) ) );
%v = @(x,c)( u(u( x, c(1:12) ),c(13:24)) );
%f = @(x,c)( cos(c(1)+v(x,c(2:25))) );
%a = @(x,c)( 1.0 + f(x,c) );
%
x = 100+linspace(-1,1,10000)*100;
f0Vals = cos( 2*pi*rf(1) + 20.0*abs(rg(1))*s3(x,2*pi*rf(7),abs(rg(7))) );
f1Vals = cos( 2*pi*rf(1) + 20.0*abs(rg(1))*s3(x,2*pi*rf(7),abs(rg(7))) );
%bigAVals = (1.0 + 0.9*cos( 2*pi*rf(1) + abs(rg(1))*s3(x,2*pi*rf(7),abs(rg(7))) )).^1;
%gVals = bigAVals.*f0Vals./(1.0+0.5*f1Vals);
gVals = f0Vals./(1.0+0.8*f1Vals);
%x = linspace(0,1,10000);
%plot( x, s(x), "o-" );
%plot( x, s1(x,rf(1),rg(1)), "o-" );
%plot( x, s1(x,2*pi*rf(1),abs(rg(1))), "o-" );
%plot( x, s2(x,rf(3),rg(3)), "o-" );
%plot( x, s2(x,2*pi*rf(3),abs(rg(3))), "o-" );
%plot( x, s2(x,rf(7),rg(7)), "o-" );
%plot( x, s3(x,2*pi*rf(7),abs(rg(7))), "o-" );
plot( x, gVals, "o-" );
grid on;
return;
