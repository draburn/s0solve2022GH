myclear;
%
sizeX = 2;
sizeF = 2;
%
if (1)
	funcPrm.matM = eye(sizeF,sizeX);
	funcPrm.vecV = [1;0];
	funcPrm.vecA = [2;0];
	funcPrm.cQuad = -1.5;
	vecXRoot = [-0.2;0];
	ax = [-0.3,1.2,-0.5,0.5];
	vecXM = [ 0.789; 0.0 ]; % From manual search.
else
	funcPrm.matM = [1,0;1,-1];
	funcPrm.vecV = [1;1];
	funcPrm.vecA = [2;1];
	funcPrm.cQuad = -1.3;
	vecXRoot = [-0.5;0.2];
	ax = [-1,1,-0.3,1.2];
	vecXM = [ 0.027; 0.683 ]; % From manual search.
end
vecFPreRoot = oneComponentCubic(vecXRoot,funcPrm);
funchF = @(x)( oneComponentCubic(x,funcPrm) - repmat(vecFPreRoot,[1,size(x,2)]) );
funchF1XY = @(x,y)( funchF([x;y])(1,:) );
funchF2XY = @(x,y)( funchF([x;y])(2,:) );
funchFNXY = @(x,y)( sqrt(sum(funchF([x;y]).^2,1)) );

if (0)
	rvecX = linspace(ax(1),ax(2),1001);
	figIndex++; figure(figIndex);
	plot( ...
	  rvecX, funchF1XY(rvecX,0*rvecX), 'o-', ...
	  rvecX, funchF2XY(rvecX,0*rvecX), 'x-' );
	grid on;
end

figIndex++; figure(figIndex);
funchVizXY = @(x,y)( asinh(funchFNXY(x,y)*10.0)/10.0 );
contourfunch( funchVizXY, ax );

if (0)
	% This space is for manual searcing of bad loc min.
	figIndex++; figure(figIndex);
	funchVizXY = @(x,y)( asinh(funchFNXY(x,y)*10.0)/10.0 );
	%contourfunch(funchVizXY,[0.788,0.789,-0.0005,0.0005]); % CASE 1
	contourfunch(funchVizXY,[0.026,0.028,0.683,0.684]); % CASE 2
return;
end

vecFM = funchF(vecXM);
vecFMHat = vecFM/norm(vecFM);
matP = eye(sizeF,sizeF) - (vecFMHat*(vecFMHat'));
funchPF = @(x)( matP * funchF(x) );
%
funchPF1XY = @(x,y)( funchPF([x;y])(1,:) );
funchPF2XY = @(x,y)( funchPF([x;y])(2,:) );
funchPFNXY = @(x,y)( sqrt(sum(funchPF([x;y]).^2,1)) );

figIndex++; figure(figIndex);
funchVizXY = @(x,y)( asinh(funchPFNXY(x,y)*10.0)/10.0 );
contourfunch( funchVizXY, ax );
