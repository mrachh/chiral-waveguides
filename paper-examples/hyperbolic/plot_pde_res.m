clear
load('everything_de1.mat');

figure(1); clf; h=pcolor(XX,YY,abs(ccur));
set(h,'EdgeColor','none');
axis equal
xlim([min(XX(:)),max(XX(:))]);
ylim([min(YY(:)),max(YY(:))]);
colorbar();
hold on;
scatter(xatoms,yatoms,40,'filled','white');
hold off
saveas(gcf,'../results/ref_hyperbolic.pdf');