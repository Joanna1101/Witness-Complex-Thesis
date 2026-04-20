% Klein Bottle


tpi = 2*pi;
npts = 100;
[theta,phi]= meshgrid(0:tpi/(npts-1):tpi,0:tpi/(npts-1):tpi);

[x,y,z,w] = KleinBottle(theta,phi);

figure();
surf(x,y,z, 'FaceAlpha',0.5,'Facecolor','interp') %'EdgeColor','none');
xlabel('x');
ylabel('y');
zlabel('z');
figure();
surf(x,z,w, 'FaceAlpha',0.5);
xlabel('x');
ylabel('z');
zlabel('w');


%Following Dłotko1, arXiv:2301.06753v1  See https://doi.org/10.1137/23M1594728
function [x,y,z,w] = KleinBottle(theta,phi)
global tpi
 x = cos(0.5.*theta).*cos(phi) - sin(0.5.*theta).*sin(2.*phi);
 y = sin(0.5.*theta).*cos(phi) + cos(0.5.*theta).*sin(2.*phi);
 z = 8.*cos(theta).*(1+ 0.5.*sin(0.5.*phi));
 w = 8.*sin(theta).*(1.0+ 0.5.*sin(0.5.*phi));
end

