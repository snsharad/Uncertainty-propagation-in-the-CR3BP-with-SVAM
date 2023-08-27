% Error v. Mahalanobis distance plot

clear;clc

load ExpM6b_case1.mat

LU = 384400;
TU = 375190;

ErrPos1 = ErrPos * LU;
MahalD1 = MahalD;
tspan11 = tspan;

clearvars -except LU TU ErrPos1 MahalD1 tspan11

load ExpM6b_case2.mat

ErrPos2 = ErrPos * LU;
MahalD2 = MahalD;
tspan2 = tspan;

clearvars -except LU TU ErrPos1 MahalD1 tspan11 ErrPos2 MahalD2 tspan2

load ExpM6b_case3.mat

ErrPos3 = ErrPos * LU;
MahalD3 = MahalD;
tspan3 = tspan;

%% Plot

Md1 = sqrt(MahalD1(1,:));
Md2 = sqrt(MahalD2(1,:));
Md3 = sqrt(MahalD3(1,:));

% id1 = 230;
% id2 = 500;
% id3 = 1000;

id1 = 10;
id2 = 300;
id3 = 900;

id11 = 11;
id12 = 310;
id13 = 900;

minP1 = min(min(ErrPos1(id11,:))); minP2 = min(min(ErrPos2(id1,:))); minP3 = min(min(ErrPos3(id1,:))); 
maxP1 = max(max(ErrPos1(id13,:))); maxP2 = max(max(ErrPos2(id3,:))); maxP3 = max(max(ErrPos3(id3,:)));

minP = min([minP1, minP2, minP3]);
maxP = max([maxP1, maxP2, maxP3]);

minColorLimit = minP;
maxColorLimit = maxP;


figure 
t = tiledlayout(3,3, 'TileSpacing', 'tight');

nexttile
% id = 50;
TF = tspan11(id11) * TU/3600;
scatter(Md1', sqrt(MahalD1(id11,:))', [], ErrPos1(id11,:), '.')
caxis([minP maxP]);
set(gca,'ColorScale','log')
axis equal
xlabel('$M_d (t_0)$', 'interpreter', 'latex')
ylabel('$M_d (t_f)$', 'interpreter', 'latex')
title("Case 1, $t_f = $" + ceil(TF) + " hr", 'interpreter', 'latex')
ylim([0 6])
xlim([0 4])
xticks(0:1:4)
yticks(0:1:6)
grid on


nexttile
% id = 250;
TF = tspan11(id12) * TU/3600;
scatter(Md1', sqrt(MahalD1(id12,:))', [], ErrPos1(id12,:), '.')
caxis([minP maxP]);
set(gca,'ColorScale','log')
axis equal
xlabel('$M_d (t_0)$', 'interpreter', 'latex')
ylabel('$M_d (t_f)$', 'interpreter', 'latex')
title("Case 1, $t_f = $" + floor(TF) + " hrs", 'interpreter', 'latex')
ylim([0 6])
xlim([0 4])
xticks(0:1:4)
yticks(0:1:6)
grid on

nexttile
% id = 500;
TF = tspan11(id13) * TU/3600;
scatter(Md1', sqrt(MahalD1(id13,:))', [], ErrPos1(id13,:), '.')
caxis([minP maxP]);
set(gca,'ColorScale','log')
axis equal
xlabel('$M_d (t_0)$', 'interpreter', 'latex')
ylabel('$M_d (t_f)$', 'interpreter', 'latex')
title("Case 1, $t_f = $" + ceil(TF) + " hrs", 'interpreter', 'latex')
ylim([0 6])
xlim([0 4])
xticks(0:1:4)
yticks(0:1:6)
grid on

%%%%%%%%%%%%%%%%%%

nexttile
% id = 50;
TF = tspan2(id1) * TU/3600;
scatter(Md2', sqrt(MahalD2(id1,:))', [], ErrPos2(id1,:), '.')
caxis([minP maxP]);
set(gca,'ColorScale','log')
axis equal
xlabel('$M_d (t_0)$', 'interpreter', 'latex')
ylabel('$M_d (t_f)$', 'interpreter', 'latex')
title({"  ";" ";"Case 2, $t_f = $" + ceil(TF) + " hr"}, 'interpreter', 'latex')
ylim([0 6])
xlim([0 4])
xticks(0:1:4)
yticks(0:1:6)
grid on


nexttile
id = 250;
TF = tspan2(id2) * TU/3600;
scatter(Md2', sqrt(MahalD2(id2,:))', [], ErrPos2(id2,:), '.')
caxis([minP maxP]);
set(gca,'ColorScale','log')
axis equal
xlabel('$M_d (t_0)$', 'interpreter', 'latex')
ylabel('$M_d (t_f)$', 'interpreter', 'latex')
title("Case 2, $t_f = $" + ceil(TF) + " hrs", 'interpreter', 'latex')
ylim([0 6])
xlim([0 4])
xticks(0:1:4)
yticks(0:1:6)
grid on

nexttile
% id = 500;
TF = tspan2(id3) * TU/3600;
scatter(Md2', sqrt(MahalD2(id3,:))', [], ErrPos2(id3,:), '.')
caxis([minP maxP]);
set(gca,'ColorScale','log')
axis equal
xlabel('$M_d (t_0)$', 'interpreter', 'latex')
ylabel('$M_d (t_f)$', 'interpreter', 'latex')
title("Case 2, $t_f = $" + ceil(TF) + " hrs", 'interpreter', 'latex')
ylim([0 6])
xlim([0 4])
xticks(0:1:4)
yticks(0:1:6)
grid on

%%%%%%%%%%%%%%%%%

nexttile
% id = 50;
TF = tspan3(id1) * TU/3600;
scatter(Md3', sqrt(MahalD3(id1,:))', [], ErrPos3(id1,:), '.')
caxis([minP maxP]);
set(gca,'ColorScale','log')
axis equal
xlabel('$M_d (t_0)$', 'interpreter', 'latex')
ylabel('$M_d (t_f)$', 'interpreter', 'latex')
title("Case 3, $t_f = $" + ceil(TF) + " hr", 'interpreter', 'latex')
ylim([0 6])
xlim([0 4])
xticks(0:1:4)
yticks(0:1:6)
grid on

nexttile
% id = 250;
TF = tspan3(id2) * TU/3600;
scatter(Md3', sqrt(MahalD3(id2,:))', [], ErrPos3(id2,:), '.')
caxis([minP maxP]);
set(gca,'ColorScale','log')
axis equal
xlabel('$M_d (t_0)$', 'interpreter', 'latex')
ylabel('$M_d (t_f)$', 'interpreter', 'latex')
title("Case 3, $t_f = $" + ceil(TF) + " hrs", 'interpreter', 'latex')
ylim([0 6])
xlim([0 4])
xticks(0:1:4)
yticks(0:1:6)
grid on


nexttile
% id = 500;
TF = tspan3(id3) * TU/3600;
scatter(Md3', sqrt(MahalD3(id3,:))', [], ErrPos3(id3,:), '.')
caxis([minP maxP]);
set(gca,'ColorScale','log')
axis equal
xlabel('$M_d (t_0)$', 'interpreter', 'latex')
ylabel('$M_d (t_f)$', 'interpreter', 'latex')
title("Case 3, $t_f = $" + ceil(TF) + " hrs", 'interpreter', 'latex')
ylim([0 6])
xlim([0 4])
xticks(0:1:4)
yticks(0:1:6)
grid on

h = colorbar;
ylabel(h, 'Error (km)')
h.Layout.Tile = 'east';
set(findall(gcf,'-property','FontSize'),'FontSize',13)


