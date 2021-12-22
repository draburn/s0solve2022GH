clear;
beta = 0.1;
gamma = linspace( 0.0, 1.0, 10000 );
phi = log(abs(gamma)./abs(beta-gamma)) ./ log(abs(beta-gamma)./abs(1.0-gamma));
plot( gamma, log(abs(phi)), 'o-' );
grid on;
axis([-0.1 1.1 -20.0 20.0]);
