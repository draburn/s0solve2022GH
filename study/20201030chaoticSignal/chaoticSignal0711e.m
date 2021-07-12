clear;
tic();
setprngstates
%setprngstates(0);
%setprngstates(78886496);
%setprngstates(20859712);
numFigs = 0;
thisFile = "chaoticSignal0711e";
%
xNormyLo = -100.0;
xNormyHi = 100.0;
rfn = @(n) rand(n,1);
rgn = @(n) randn(n,1);
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
f = @(x,rf,rg) cos( 2*pi*rf(1) + 20.0*abs(rg(1))*s3( x, 2*pi*rf(2:8), abs(rg(2:8)) ) );
%
x = 100+linspace(-1,1,10000)*100;
%f0Vals = cos( 2*pi*rfn(1) + 20.0*abs(rgn(1))*s3(x,2*pi*rfn(7),abs(rgn(7))) );
%f1Vals = cos( 2*pi*rfn(1) + 20.0*abs(rgn(1))*s3(x,2*pi*rfn(7),abs(rgn(7))) );
%gVals = f0Vals./(1.0+0.8*f1Vals);
plot( x, f(x,rfn(8),rgn(8))./(1.0+0.8*f(x,rfn(8),rgn(8))), "o-" );
grid on;
%
toc;
thisFile = [ "RETURN FROM " thisFile ];
return;
