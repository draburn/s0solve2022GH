numFigs = 0;
funchVizLog = @(x)(log( eps025*max(max(x)) + x - min(min(x)) ));
%
%
%
numFigs++; figure(numFigs);
contourf( x1Mesh, x2Mesh, funchVizLog(omegaMesh), 30 );
hold on;
plot( matY_grad_old(1,:), matY_grad_old(2,:), 'rx-', 'linewidth', 2, 'markersize', 5 );
plot( matY_grad_simple(1,:), matY_grad_simple(2,:), 'go-', 'linewidth', 2, 'markersize', 5 );
plot( matY_grad1(1,:), matY_grad1(2,:), 'b+-', 'linewidth', 2, 'markersize', 5 );
plot( matY_ocqLev(1,:), matY_ocqLev(2,:), 'y*-', 'linewidth', 2, 'markersize', 5 );
hold off;
grid on;
if (useAxisEqual)
	axis equal;
end
colormap( mycmap );
xlabel( "x1" );
ylabel( "x2" );
title( "log(omega-omegaMin) vs x1, x2" );
%
%
%
numFigs++; figure(numFigs);
semilogy( ...
  sqrt(sum(diff(matY_grad_old').^2,2)), 'rx-', ...
  sqrt(sum(diff(matY_grad_simple').^2,2)), 'go-', ...
  sqrt(sum(diff(matY_grad1').^2,2)), 'b+-', ...
  sqrt(sum(diff(matY_ocqLev').^2,2)), 'y*-' );
grid on;
xlabel( "step index" );
ylabel( "step size" );
title( "step size vs step index" );
%
%
%
numFigs++; figure(numFigs);
semilogy( ...
  0.5*sqrt(sum((matF_grad_old').^2,2)), 'rx-', ...
  0.5*sqrt(sum((matF_grad_simple').^2,2)), 'go-', ...
  0.5*sqrt(sum((matF_grad1').^2,2)), 'b+-', ...
  0.5*sqrt(sum((matF_ocqLev').^2,2)), 'y*-' );
grid on;
xlabel( "step index" );
ylabel( "omega" );
title( "omega vs step index" );
%
%
%
numFigs++; figure(numFigs);
plot( ...
  sqrt(sum(matD_grad_old.^2,1)), 'rx-', ...
  sqrt(sum(matD_grad_simple.^2,1)), 'go-', ...
  sqrt(sum(matD_grad1.^2,1)), 'b+-', ...
  sqrt(sum(matD_ocqLev.^2,1)), 'y*-' );
grid on;
xlabel( "step index" );
ylabel( "distance" );
title( "distance vs step index" );
%
%
%
numFigs++; figure(numFigs);
semilogy( ...
  sqrt(sum(matD_grad_old.^2,1)), 0.5*sqrt(sum(matF_grad_old.^2,1)), 'rx-', ...
  sqrt(sum(matD_grad_simple.^2,1)), 0.5*sqrt(sum(matF_grad_simple.^2,1)), 'go-', ...
  sqrt(sum(matD_grad1.^2,1)), 0.5*sqrt(sum(matF_grad1.^2,1)), 'b+-', ...
  sqrt(sum(matD_ocqLev.^2,1)), 0.5*sqrt(sum(matF_ocqLev.^2,1)), 'y*-' );
grid on;
xlabel( "distance" );
ylabel( "omega" );
title( "omega vs distance" );
%
%
