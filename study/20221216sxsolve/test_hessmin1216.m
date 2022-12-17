clear;
setprngstates(0);
sizeX = 5

fRoot = 0.0;
matH = mtm(randn(sizeX-1,sizeX))
vecXRoot = randn(sizeX,1)
funchD = @(x)( x - vecXRoot );
funchHDOfD = @(d)( matH*d );
funchFOfD = @(d)( fRoot + 0.5*sum( d .* funchHDOfD(d), 1 ) );
funchGOfD = @(d)( funchHDOfD(d) );
funchF = @(x)( funchFOfD(funchD(x)) );
funchG = @(x)( funchGOfD(funchD(x)) );

vecX0 = zeros(sizeX,1);
f0 = funchF(vecX0)
vecG0 = funchG(vecX0)

[ vecDeltaCalc, dat ] = hessmin1216( f0, vecG0, matH )
fNext = funchF(vecX0+vecDeltaCalc)