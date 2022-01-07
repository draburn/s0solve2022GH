% Sample usage...
%   [ gridX1, gridX2, gridF ] = genVizGrids( funch F );
%   function [ gridX1, gridX2, gridF, gridCX1, gridCX2, gridD1F, gridD2F, gridMagDF, gridThetaDF ] = genVizGrids( ...
%     funchF, ...
%     isVectorized = false, ...
%     ax = [ -3.0, 3.0, -5.0, 5.0 ], ...
%     numXVals = [ 22, 26 ], ...
%     prm = [] )

function [ gridX1, gridX2, gridF, gridCX1, gridCX2, gridD1F, gridD2F, gridMagDF, gridThetaDF ] = genVizGrids( ...
  funchF, ...
  isVectorized = false, ...
  ax = [ -5.0, 5.0, -5.0, 5.0 ], ...
  numXVals = [ 22, 26 ], ...
  prm = [] )
	%
	x1UniqVals = linspace( ax(1), ax(2), numXVals(1) );
	x2UniqVals = linspace( ax(3), ax(4), numXVals(2) );
	[ gridX1, gridX2 ] = ndgrid( x1UniqVals, x2UniqVals );
	%
	gridF = zeros( numXVals(1), numXVals(2) );
	if ( ~isVectorized )
		parfor i1 = 1:numXVals(1)
		parfor i2 = 1:numXVals(2)
			gridF(i1,i2) = funchF( [x1UniqVals(i1);x2UniqVals(i2)] );
		end
		end
	else
		vecXVals = [ reshape( gridX1, 1, [] ); reshape( gridX2, 1, [] ) ];
		vecFVals = funchF( vecXVals );
		gridF = reshape(vecFVals,[numXVals(1),numXVals(2)]);
		clear vecXVals;
		clear vecFVals;
	end
	%
	if ( nargout >= 4 )
		gridCX1 = gridX1(2:end-1,2:end-1);
		gridCX2 = gridX2(2:end-1,2:end-1);
		gridD1F = (gridF(3:end,2:end-1)-gridF(1:end-2,2:end-1))./(x1UniqVals(3:end)-x1UniqVals(1:end-2))';
		gridD2F = (gridF(2:end-1,3:end)-gridF(2:end-1,1:end-2))./(x2UniqVals(3:end)-x2UniqVals(1:end-2));
		if ( nargout >= 8 )
			gridMagDF = sqrt( gridD1F.^2 + gridD2F.^2 );
			gridThetaDF = atan2( gridD2F, gridD1F );
		end
	end
return;
end


%!test
%!	thisFile = "genVizGrids test";
%!	setprngstates();
%!	numFigs = 0;
%!	%
%!	sizeX = 2;
%!	switch (2)
%!	case 1
%!		sizeF = sizeX;
%!		vecXCent = zeros( sizeX, 1 );
%!		matA = eye( sizeF, sizeX );
%!	case 2
%!		sizeF = 2;
%!		vecXCent = randn( sizeX, 1 );
%!		matA = eye(sizeF,sizeX);%randn( sizeF, sizeX );
%!	end
%!	funchF = @(x)( 0.5*sumsq(matA*(x-vecXCent)) );
%!	%
%!	isVectorized = false;
%!	[ gridX1, gridX2, gridF, gridCX1, gridCX2, gridD1F, gridD2F, gridMagDF, gridThetaDF ] = genVizGrids( funchF, isVectorized ); 
%!	%[ gridX1, gridX2, gridF ] = genVizGrids( funchF, isVectorized );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridX1, gridX2, sqrt(gridF) );
%!	grid on;
%!	title( "sqrt(F) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	isVectorized = true;
%!	[ gridX1, gridX2, gridF, gridCX1, gridCX2, gridD1F, gridD2F, gridMagDF, gridThetaDF ] = genVizGrids( funchF, isVectorized ); 
%!	%[ gridX1, gridX2, gridF ] = genVizGrids( funchF, isVectorized );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridX1, gridX2, sqrt(gridF) );
%!	grid on;
%!	title( "sqrt(F) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridCX1, gridCX2, gridD1F );
%!	grid on;
%!	title( "dF/dx1 vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridCX1, gridCX2, gridD2F );
%!	grid on;
%!	title( "dF/dx2 vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridCX1, gridCX2, sqrt(abs(gridD1F)) );
%!	grid on;
%!	title( "sqrt(|dF/dx1|) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridCX1, gridCX2, sqrt(abs(gridD2F)) );
%!	grid on;
%!	title( "sqrt(|dF/dx2|) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridCX1, gridCX2, sqrt(gridMagDF) );
%!	grid on;
%!	title( "||nablaF|| vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	contourf( gridCX1, gridCX2, gridThetaDF );
%!	grid on;
%!	title( "theta(nablaF) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	numFigs++; figure(numFigs);
%!	%contourf( gridCX1, gridCX2, gridThetaDF );
%!	imagesc( gridCX1', gridCX2', gridThetaDF' );
%!	set(gca,'ydir','normal');
%!	colormap(hsv);
%!	grid on;
%!	title( "theta(nablaF) vs (x1,x2)" );
%!	xlabel( "x1" );
%!	ylabel( "x2" );
%!	%
%!	msg( thisFile, __LINE__, "*** Please manually confirm the figures look correct. ***" );
%!	%
%!	return
%!	%
%!	size(gridCX1)
%!	gridCX1( 8:12, 10:14 )
%!	size(gridCX2)
%!	gridCX2( 8:12, 10:14 )
%!	%
%!	size(gridD1F)
%!	gridD1F( 8:12, 10:14 )
%!	size(gridD2F)
%!	gridD2F( 8:12, 10:14 )
%!	%
%!	size(gridThetaDF)
%!	gridThetaDF( 8:12, 10:14 )
