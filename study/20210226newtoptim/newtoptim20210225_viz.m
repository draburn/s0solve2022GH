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
toc();
