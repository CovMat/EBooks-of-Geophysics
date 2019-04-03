% plot the "L-curve" : different Damp with a certain Smooth
% Hao Guo. 2016.7.8

close all
clc
clear

%% load solution norm, residual norm with different Parameters
curve = load('lcurve.dat');


%%  Normalized solution (slowness) norm versus Normalizedabsolute residual norm
figure
range1=max(curve(:,3))-min(curve(:,3));
range2=max(curve(:,5))-min(curve(:,5));
plot((curve(:,3)-min(curve(:,3)))/range1,(curve(:,5)-min(curve(:,5)))/range2,'r-*');

text((curve(:,3)-min(curve(:,3)))/range1-0.00,(curve(:,5)-min(curve(:,5)))/range2+0.00,num2str(curve(:,2)));

xlabel('Normalized solution norm');ylabel('Normalized absolute residual norm');
axis equal;
axis([0,1,0,1]);
print(gcf,'-djpeg','normalizedlcurve_ChooseDamp_absres','-r720');

%%  Normalized solution (slowness) norm versus Normalized weighted residual norm
figure
range1=max(curve(:,3))-min(curve(:,3));
range2=max(curve(:,6))-min(curve(:,6));
plot((curve(:,3)-min(curve(:,3)))/range1,(curve(:,6)-min(curve(:,6)))/range2,'r-*');

text((curve(:,3)-min(curve(:,3)))/range1-0.00,(curve(:,6)-min(curve(:,6)))/range2+0.0,num2str(curve(:,2)));

xlabel('Normalized solution norm');ylabel('Normalized weighted residual norm');
axis equal;
axis([0,1,0,1]);

print(gcf,'-djpeg','normalizedlcurve_ChooseDamp_wtres','-r720');

%% Unnormalized solution (slowness) norm versus Unnormalized absolute residual norm
figure
plot(curve(:,3),curve(:,5),'r-*');

text(curve(:,3),curve(:,5),num2str(curve(:,2)));
xlabel('solution norm');ylabel('absolute residual norm');
print(gcf,'-djpeg','lcurve_ChooseDamp_absres','-r720');

%% Unnormalized solution (slowness) norm versus Unnormalized weighted residual norm
figure
plot(curve(:,3),curve(:,6),'r-*');

text(curve(:,3),curve(:,6),num2str(curve(:,2)));
xlabel('solution norm');ylabel('weighted residual norm');
print(gcf,'-djpeg','lcurve_ChooseDamp_wtres','-r720');
%