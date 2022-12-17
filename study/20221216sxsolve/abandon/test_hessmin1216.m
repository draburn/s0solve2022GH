printf( "\n\n\nvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n" );
msg( __FILE__, __LINE__, "Hai!" );
clear;
setprngstates(0);
sizeX = 5

fRoot = 0.0;
vecXRoot = randn(sizeX,1);
vecX0 = zeros(sizeX,1);
vecD0 = vecX0-vecXRoot;
[ norm(vecD0), 0.0 ]
%
msg( __FILE__, __LINE__, "Strongly pos-def..." );
matH = mtm(randn(sizeX,sizeX));
f0 = fRoot + 0.5*vecD0'*matH*vecD0;
vecG0 = matH*vecD0;
[ vecDeltaCalc, dat ] = hessmin1216( f0, vecG0, matH );
vecXNext = vecX0 + vecDeltaCalc;
vecDNext = vecXNext - vecXRoot;
[ norm(vecDNext), fRoot + 0.5*vecDNext'*matH*vecDNext ]
%
msg( __FILE__, __LINE__, "Weakly pos-def..." );
matH = mtm(randn(1,sizeX)) + sqrt(eps)*mtm(randn(sizeX,sizeX));
f0 = fRoot + 0.5*vecD0'*matH*vecD0;
vecG0 = matH*vecD0;
[ vecDeltaCalc, dat ] = hessmin1216( f0, vecG0, matH );
vecXNext = vecX0 + vecDeltaCalc;
vecDNext = vecXNext - vecXRoot;
[ norm(vecDNext), fRoot + 0.5*vecDNext'*matH*vecDNext ]
%
msg( __FILE__, __LINE__, "Weakly pos-semi-def..." );
matH = mtm(randn(sizeX-1,sizeX));
f0 = fRoot + 0.5*vecD0'*matH*vecD0;
vecG0 = matH*vecD0;
[ vecDeltaCalc, dat ] = hessmin1216( f0, vecG0, matH );
vecXNext = vecX0 + vecDeltaCalc;
vecDNext = vecXNext - vecXRoot;
[ norm(vecDNext), fRoot + 0.5*vecDNext'*matH*vecDNext ]
%
msg( __FILE__, __LINE__, "Strongly pos-semi-def..." );
matH = mtm(randn(1,sizeX));
f0 = fRoot + 0.5*vecD0'*matH*vecD0;
vecG0 = matH*vecD0;
[ vecDeltaCalc, dat ] = hessmin1216( f0, vecG0, matH );
vecXNext = vecX0 + vecDeltaCalc;
vecDNext = vecXNext - vecXRoot;
[ norm(vecDNext), fRoot + 0.5*vecDNext'*matH*vecDNext ]
%
fRoot = 100.0
%
msg( __FILE__, __LINE__, "Weakly has-neg..." );
matH = mtm(randn(sizeX-1,sizeX)) - sqrt(eps)*mtm(randn(sizeX,sizeX));
f0 = fRoot + 0.5*vecD0'*matH*vecD0;
vecG0 = matH*vecD0;
[ vecDeltaCalc, dat ] = hessmin1216( f0, vecG0, matH );
vecXNext = vecX0 + vecDeltaCalc;
vecDNext = vecXNext - vecXRoot;
[ norm(vecDNext), fRoot + 0.5*vecDNext'*matH*vecDNext ]

msg( __FILE__, __LINE__, "Strongly has-neg..." );
matH = mtm(randn(sizeX-1,sizeX)) - mtm(randn(1,sizeX));
f0 = fRoot + 0.5*vecD0'*matH*vecD0;
vecG0 = matH*vecD0;
[ vecDeltaCalc, dat ] = hessmin1216( f0, vecG0, matH );
vecXNext = vecX0 + vecDeltaCalc;
vecDNext = vecXNext - vecXRoot;
[ norm(vecDNext), fRoot + 0.5*vecDNext'*matH*vecDNext ]
%
msg( __FILE__, __LINE__, "Fully neg-def..." );
matH = -mtm(randn(sizeX,sizeX));
f0 = fRoot + 0.5*vecD0'*matH*vecD0;
vecG0 = matH*vecD0;
[ vecDeltaCalc, dat ] = hessmin1216( f0, vecG0, matH );
vecXNext = vecX0 + vecDeltaCalc;
vecDNext = vecXNext - vecXRoot;
[ norm(vecDNext), fRoot + 0.5*vecDNext'*matH*vecDNext ]
