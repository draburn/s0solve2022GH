tic();
msg("newtoptim20210225_viz",__LINE__,"Generating visualizations...");
numFigs = 0;
%
numFigs++; figure(numFigs);
plot( vecX, vecProbLETau(:,1), 'o-' );
hold on;
for n=2:length(list_taus)
	plot( vecX, vecProbLETau(:,n), 'o-' );
end
plot( vecX, 0*vecX, 'k-' );
hold off;
grid on;
xlabel( "x" );
ylabel( "P[|F|<=tau]" );
title( "P[|F|<=tau] vs x" );
%
numFigs++; figure(numFigs);
plot( vecX, vecF_avgSignPows(:,1), 'o-' );
hold on;
for n=2:length(list_pows)
	plot( vecX, vecF_avgSignPows(:,n), 'o-' );
end
plot( vecX, 0*vecX, 'k-' );
hold off;
grid on;
xlabel( "x" );
ylabel( "<F^p>^1^/^p" );
title( "<F^p>^1^/^p vs x" );
%
numFigs++; figure(numFigs);
plot( vecX, vecF_avgAbsPows(:,1), 'o-' );
hold on;
for n=2:length(list_pows)
	plot( vecX, vecF_avgAbsPows(:,n), 'o-' );
end
plot( vecX, 0*vecX, 'k-' );
hold off;
grid on;
xlabel( "x" );
ylabel( "<|F|^p>^1^/^p" );
title( "<|F|^p>^1^/^p vs x" );
%
numFigs++; figure(numFigs);
plot( vecX, vecF_signPercentiles(:,1), 'o-' );
hold on;
for n=2:length(list_percentiles)
	plot( vecX, vecF_signPercentiles(:,n), 'o-' );
end
plot( vecX, 0*vecX, 'k-' );
hold off;
grid on;
xlabel( "x" );
ylabel( "F percentiles" );
title( "F percentiles vs x" );
%
numFigs++; figure(numFigs);
plot( vecX, vecF_absPercentiles(:,1), 'o-' );
hold on;
for n=2:length(list_percentiles)
	plot( vecX, vecF_absPercentiles(:,n), 'o-' );
end
plot( vecX, 0*vecX, 'k-' );
hold off;
grid on;
xlabel( "x" );
ylabel( "|F| percentiles" );
title( "|F| percentiles vs x" );
%
%
% d/dx...
numFigs++; figure(numFigs);
plot( cent(vecX), diff(vecProbLETau(:,1))./diff(vecX), 'o-' );
hold on;
for n=2:length(list_taus)
	plot( cent(vecX), diff(vecProbLETau(:,n))./diff(vecX), 'o-' );
end
plot( vecX, 0*vecX, 'k-' );
hold off;
grid on;
xlabel( "x" );
ylabel( "d/dx P[|F|<=tau]" );
title( "d/dx P[|F|<=tau] vs x" );
%
numFigs++; figure(numFigs);
plot( cent(vecX), diff(vecF_avgSignPows(:,1))./diff(vecX), 'o-' );
hold on;
for n=2:length(list_pows)
	plot( cent(vecX), diff(vecF_avgSignPows(:,n))./diff(vecX), 'o-' );
end
plot( vecX, 0*vecX, 'k-' );
hold off;
grid on;
xlabel( "x" );
ylabel( "d/dx <F^p>^1^/^p" );
title( "d/dx <F^p>^1^/^p vs x" );
%
numFigs++; figure(numFigs);
plot( cent(vecX), diff(vecF_avgAbsPows(:,1))./diff(vecX), 'o-' );
hold on;
for n=2:length(list_pows)
	plot( cent(vecX), diff(vecF_avgAbsPows(:,n))./diff(vecX), 'o-' );
end
plot( vecX, 0*vecX, 'k-' );
hold off;
grid on;
xlabel( "x" );
ylabel( "d/dx <|F|^p>^1^/^p" );
title( "d/dx <|F|^p>^1^/^p vs x" );
%
numFigs++; figure(numFigs);
plot( cent(vecX), diff(vecF_signPercentiles(:,1))./diff(vecX), 'o-' );
hold on;
for n=2:length(list_percentiles)
	plot( cent(vecX), diff(vecF_signPercentiles(:,n))./diff(vecX), 'o-' );
end
plot( vecX, 0*vecX, 'k-' );
hold off;
grid on;
xlabel( "x" );
ylabel( "d/dx F percentiles" );
title( "d/dx F percentiles vs x" );
%
numFigs++; figure(numFigs);
plot( cent(vecX), diff(vecF_absPercentiles(:,1))./diff(vecX), 'o-' );
hold on;
for n=2:length(list_percentiles)
	plot( cent(vecX), diff(vecF_absPercentiles(:,n))./diff(vecX), 'o-' );
end
plot( vecX, 0*vecX, 'k-' );
hold off;
grid on;
xlabel( "x" );
ylabel( "d/dx |F| percentiles" );
title( "d/dx |F| percentiles vs x" );
%
%
for n=1:length(list_atXs)
	numFigs++; figure(numFigs);
	hist( fVals_signAtX(n,:), linspace(-4,1,51) );
	axis([-2,1,0,150000]);
	grid on;
	title( ["Histogram for F at index " num2str(n)] );
end
%
for n=1:length(list_atXs)
	numFigs++; figure(numFigs);
	hist( fVals_absAtX(n,:), linspace(-0.01,2,101) );
	axis([0,1.5,0,40000]);
	grid on;
	title( ["Histogram for |F| at index " num2str(n)] );
end
%
%
cumlhist_numPtsDesired = 1000;
cumlhist_stride = round(numTrials/cumlhist_numPtsDesired);
cumlhist_numPts = floor(numTrials/cumlhist_stride);
tVals = linspace(0.0,1.0,cumlhist_numPts);
numFigs++; figure(numFigs);
plot( fVals_signAtX(1,1:cumlhist_stride:end), tVals, 'o-' );
hold on;
for n=2:length(list_atXs)
	plot( fVals_signAtX(n,1:cumlhist_stride:end), tVals, 'o-' );
end
hold off;
axis([-2.5,1.5,0,1.1]);
grid on;
%
numFigs++; figure(numFigs);
plot( fVals_absAtX(1,1:cumlhist_stride:end), tVals, 'o-' );
hold on;
for n=2:length(list_atXs)
	plot( fVals_absAtX(n,1:cumlhist_stride:end), tVals, 'gx-' );
end
hold off;
axis([0,10,0,1]);
grid on;
%
%
toc();
