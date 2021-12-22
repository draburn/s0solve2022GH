clear;
commondefs;
thisFile = "lookatLoggishRes4_quadSym";
numFigs = 100;
%
caseNum = 1;
switch (caseNum)
case 1
	setprngstates(44236336);
	c = randn(20,1);
	funchF_pre = @(x)( c(1) + c(2)*x + 0*c(3)*x.^2 + 0*c(4)*x.^3 ...
	 + c(5) * cos( c(6) + c(7)*x ) + c(8) * cos( c(8) + c(9) * x ) ...
	 + c(10) * cos( c(11) + c(12)*cos( c(13) + c(14)*x) ) );
	xSecret = randn*exp(abs(randn));
	fPreSecret = funchF_pre(xSecret);
	funchF = @(x)( funchF_pre(x) - fPreSecret ).^2;
	x_secret = -0.890874387338734; % This is the left abs min.
	%%xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843, 0.0 ]);
	xVals = sort([ -1.00057993392590, -1.0, -0.971949673093843 ]);
	xVals = sort([ -0.831698038563058, xVals ]);
	xVals = sort([ -0.897620231966365, xVals ]);
	xVals = sort([ -0.867992251154230, xVals ]);
	xVals = sort([ -0.88284017185230, xVals ]);
	xVals = sort([ -0.887137697860992, xVals ]); % xB goes past data.
	xVals = sort([ -0.886078607798801, xVals ]);
	xVals = sort([ -0.890138678468405, xVals ]); % xB goes way past data!
	xVals = sort([ -0.889487776317036, xVals ]); % must balance...
	xVals = sort([ -0.893393189225250, xVals ]);
	xVals = sort([ -0.890922389166615, xVals ]); % xB goes past...
	xVals = sort([ -0.890775831691018, xVals ]); % must balance...
	xVals = sort([ -0.891655176544600, xVals ]);
	xVals = sort([ -0.890874186962029, xVals ]); %1E-16.
case 3
	x_secret = 0.8;
	funchF = @(x)( abs(x-x_secret).^6.5 );
	xVals = x_secret+linspace(-1.0,1.0,8)-0.3;
otherwise
	error(["Invalid value of caseNum (", num2str(caseNum), ")."]);
end
%
%
%
numPts = size(xVals,2);
fVals = funchF(xVals);
%
viz_numPts = 5000;
viz_xVals = linspace(min(xVals),max(xVals),viz_numPts);
%
%
%
xNew = quadSymInterp( xVals, fVals )
xNewRes = xNew-x_secret
%
numFigs++; figure(numFigs);
plot( ...
  xVals, fVals, 'ko', 'linewidth', 2, 'markersize', 25, ...
  viz_xVals, funchF(viz_xVals),  'o-', "color", [0.4,0.4,0.4], "markersize", 4, ...
  x_secret, funchF(x_secret), '*', "color", [0.8,0.7,0.0], "linewidth", 4, "markersize", 25, ...
  xNew, funchF(xNew), '+', "color", [0.8,0.0,1.0], "linewidth", 4, "markersize", 25, ...
  viz_xVals, 0*viz_xVals, 'k-' );
grid on;
xlabel( "x" );
ylabel( "f" );
title( "f vs x" );
