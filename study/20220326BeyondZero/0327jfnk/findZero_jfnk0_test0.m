clear;
setVerbLevs;
%setprngstates(39901520); % 39901520 oldCrude is particularly bad relative to fsovle.
setprngstates(0);
numFigs = 0;
%
sizeX = 100;
sizeF = 100;
%
vecXE = randn(sizeX,1);
matJE = diag(1.0+abs(randn(min([sizeF,sizeX]),1))) + 0.00*randn(sizeF,sizeX);
matA0 = 0.0001*randn(sizeF,sizeX);
matA1 = randn(sizeX,sizeX);
matA2 = randn(sizeX,sizeX);
matB0 = 0.0001*randn(sizeF,sizeX);
matB1 = randn(sizeX,sizeX);
matB2 = randn(sizeX,sizeX);
matB3 = randn(sizeX,sizeX);
y = @(x)( x - vecXE );
funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
msg( __FILE__, __LINE__, sprintf( "rcond(matJE'*matJE) = %f.", rcond(matJE'*matJE) ) );
%
vecX0 = zeros(sizeX,1);
vecF0 = funchF(vecX0);
%
timeSS = time();
[ vecXF_fsolve, vecFF_fsolve, datOut_fsolve ] = findZero_fsolve( vecX0, funchF );
time_fsolve = time()-timeSS;
%
timeSS = time();
[ vecXF_oldCrude, vecFF_oldCrude, datOut_oldCrude ] = findZero_oldCrude( vecX0, funchF );
time_oldCrude = time()-timeSS;
%
timeSS = time();
[ vecXF_shouldIterMin, vecFF_shouldIterMin, datOut_shouldIterMin ] = findZero_shouldIterMin( vecX0, funchF );
time_shouldIterMin = time()-timeSS;
%
msg( __FILE__, __LINE__, sprintf( "norm(vecF0) = %g.", norm(vecF0) ) );
msg( __FILE__, __LINE__, "Solver results..." );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11s;  %11s;  %11s;  %11s.", "solver name", "||vecFF||", "iterCount", "fevalCount", "time(s)" ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "fsolve", norm(vecFF_fsolve), datOut_fsolve.iterCount, datOut_fsolve.fevalCount, time_fsolve ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "oldCrude", norm(vecFF_oldCrude), datOut_oldCrude.iterCount, datOut_oldCrude.fevalCount, time_oldCrude ) );
msg( __FILE__, __LINE__, sprintf( "  %15s:  %11.3e;  %11d;  %11d;  %11.3e.", "shouldIterMin", norm(vecFF_shouldIterMin), datOut_shouldIterMin.iterCount, datOut_shouldIterMin.fevalCount, time_shouldIterMin ) );
