clear;
tic();
%
setprngstates(0);
sizeX = 5;
%
sizeF = sizeX;
vecXE = randn(sizeX,1);
matJE = diag(1.0+abs(randn(min([sizeF,sizeX]),1))) + 0.3*randn(sizeF,sizeX);
matA0 = 0.0001*randn(sizeF,sizeX);
matA1 = randn(sizeX,sizeX);
matA2 = randn(sizeX,sizeX);
matB0 = 0.0001*randn(sizeF,sizeX);
matB1 = randn(sizeX,sizeX);
matB2 = randn(sizeX,sizeX);
matB3 = randn(sizeX,sizeX);
%
y = @(x)( x - vecXE );
funchF = @(x)( matJE*y(x) + matA0*( (matA1*y(x)) .* (matA2*y(x)) ) + matB0*( (matB1*y(x)) .* (matB2*y(x)) .* (matB3*y(x)) ) );
%
vecX0 = zeros(sizeX,1);
vecF0 = funchF(vecX0);
prm = [];
msg( __FILE__, __LINE__, "Calling zlinsolf()..." );
[ vecXF, vecFF, datOut ] = zlinsolf100( funchF, vecX0, vecF0, prm );
msg( __FILE__, __LINE__, "Returned from zlinsolf()." );
vecFF = funchF(vecXF);
[ norm(vecF0), norm(vecFF) ]
%[ vecX0, vecXF ]
%[ vecF0, vecFF ]
msg( __FILE__, __LINE__, "Goodbye." );
%
toc();
