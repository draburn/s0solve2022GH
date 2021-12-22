%
%ax = [ -0.1, 0.1, -0.01, 0.01 ];

numX1Vals = 101;
numX2Vals = 101;
x1Vals = linspace(ax(1),ax(2),numX1Vals);
x2Vals = linspace(ax(3),ax(4),numX2Vals);
[ x1Mesh, x2Mesh ] = meshgrid( x1Vals, x2Vals );
matX = [ reshape(x1Mesh,1,[]); reshape(x2Mesh,1,[]) ];
matF = funchF(matX,testFuncPrm);
f1Mesh = reshape(matF(1,:),numX2Vals,numX1Vals);
f2Mesh = reshape(matF(2,:),numX2Vals,numX1Vals);
omegaMesh = 0.5*( f1Mesh.^2 + f2Mesh.^2 );
%
msg( thisFile, __LINE__, sprintf( "F1 scale: %g to %g.", min(min((f1Mesh))), max(max((f1Mesh))) ) );
msg( thisFile, __LINE__, sprintf( "F2 scale: %g to %g.", min(min((f2Mesh))), max(max((f2Mesh))) ) );
msg( thisFile, __LINE__, sprintf( "omega scale: %g to %g.", min(min((omegaMesh))), max(max((omegaMesh))) ) );
%
funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
