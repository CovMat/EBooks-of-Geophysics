close all;
clc
clear;
%
% plot travel time versus hypocenter distance.
% Input file: t_dist.dat
%
% Please determine b1_P, b2_P, b1_S, b2_S velues according to the P-wave
% and S-wave time-dist curve; 
% Please record the output information: slope_P, b_P; slope_S, b_S
% These are the parameters for pha_t-dist_selection.awk
%
% Hao Guo. 2016.7.9
%


%% Parmeters
b1_P = 7;
b2_P = 7;
b1_S = 7;
b2_S = 7;

%% separate P and S t-dist data.
t_d = load('t_dist.dat');
n1 = 0; n2 = 0;
for i = 1:length(t_d(:,1))
    if t_d(i,3)==1
        n1 = n1+1;
        t_d_P(n1,1) = t_d(i,1);
        t_d_P(n1,2) = t_d(i,2);
    elseif t_d(i,3)==2
        n2 = n2+1;
        t_d_S(n2,1) = t_d(i,1);
        t_d_S(n2,2) = t_d(i,2);
    end
end
%% plot t-dist curve of P wave
figure;
subplot(1,2,1);
plot(t_d_P(:,2),t_d_P(:,1),'r.');hold on;
p = polyfit(t_d_P(:,2),t_d_P(:,1),1);
x1 = linspace(min(t_d_P(:,2)),max(t_d_P(:,2)));
y1 = polyval(p,x1);
plot(x1,y1,'k');
p2 = p; p2(1,2) = p2(1,2) + b1_P;
p3 = p; p3(1,2) = p3(1,2) - b2_P;
y2 = polyval(p2,x1);
plot(x1,y2,'g'); hold on;
y3 = polyval(p3,x1);
plot(x1,y3,'g'); hold on;

axis([0 800 0 200]);
title('t-dist curve of P');
xlabel('Hypocenter Distance (km)');ylabel('Travel Time (s)');

fprintf('P wave: the slope and cutoff value of the fitted curve (slope_P, b_P) are %f and %f\n',p(1),p(2));
%% plot t-dist curve of S wave
subplot(1,2,2); 
plot(t_d_S(:,2),t_d_S(:,1),'r.');hold on;
p = polyfit(t_d_S(:,2),t_d_S(:,1),1);
x1 = linspace(min(t_d_S(:,2)),max(t_d_S(:,2)));
y1 = polyval(p,x1);
plot(x1,y1,'k');
p2 = p; p2(1,2) = p2(1,2) + b1_P;
p3 = p; p3(1,2) = p3(1,2) - b2_P;
y2 = polyval(p2,x1);
plot(x1,y2,'g'); hold on;
y3 = polyval(p3,x1);
plot(x1,y3,'g'); hold on;

axis([0 800 0 200]);
title('t-dist curve of S');
xlabel('Hypocenter Distance (km)');
print('-depsc2','t_dist','-r300');

fprintf('S wave: the slope and cutoff value of the fitted curve (slope_S, b_S) are %f and %f\n',p(1),p(2));