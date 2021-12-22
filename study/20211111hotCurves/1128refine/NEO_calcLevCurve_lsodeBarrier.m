function matX = NEO_calcLevCurve_lsodeBarrier( funchG, vecX0, prm=[] )
	thisFile = "NEO_calcLevCurve_lsodeBarrier";
	%
	numDVals = 51;
	dHi = 4.0;
	dVals = linspace( 0.0, dHi, numDVals );
	%barrierThickness = dHi/(5.0*(numDVals-1.0));
	barrierThickness = 0.01;
	%
	matX(:,1) = vecX0;
	for n=2:numDVals
		d = dVals(n);
		funchZ0 = @(x)( (norm(x-vecX0)-d)/barrierThickness );
		funchZ1 = @(x)( funchZ0(x).*(funchZ0(x)>0.0) );
		funchL = @(x)( exp( funchZ1(x).^4 ) - 1.0 );
		%%%funchL = @(x)( exp( funchZ1(x).^2 ) - 1.0 );
		funchAntiG = @(x,t)( -funchG(x,t)-(x-vecX0)*funchL(x) );
		matX_temp = lsode( funchAntiG, matX(:,n-1), [0.0,4.0] )';
		matX(:,n) = matX_temp(:,2);
		%[ norm(matX(:,n)-vecX0), norm(matX(:,n)-vecX0)-d, funchZ0(matX(:,n)), funchZ1(matX(:,n)), funchL(matX(:,n)) ]
	end
	%funchAntiG = @(x)( -funchG(x) );
	%matX = lsode( funchAntiG, vecX0, 1000.0*linspace(0.0,1.0,1001) )';
	%
return;
end
