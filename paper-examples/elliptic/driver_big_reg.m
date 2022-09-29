clear


A = load('../hyperbolic/data/atoms.mat');
run ../../startup.m
addpath ../../../fmm2d/matlab

dlam = A.s0;
dk   = A.de;


nover = 10;
wavelam = 2*pi/dk;

nx0 = 100;
ny0 = 100;
nx = nx0*nover;
ny = ny0*nover;

halfwave = wavelam/2/nover;
xx = -(nx-1)/2*halfwave + (0:halfwave:(nx-1)*halfwave);
yy = -(ny-1)/2*halfwave + (0:halfwave:(ny-1)*halfwave);

[xx_use,yy_use]= meshgrid(xx,yy);

xmin = min(xx);
xmax = max(xx);
ymin = min(yy);
ymax = max(yy);

dx = xmax-xmin;
dy = ymax-ymin;



xmin = xmin - dx/4;
xmax = xmax + dx/4;

ymin = ymin - dy/4;
ymax = ymax + dy/4;

xs = xmin:(halfwave/10*nover):xmax;
ys = ymin:(halfwave/10*nover):ymax;


xa = xx_use(:);
ya = yy_use(:);
natoms = length(xa);

ra = [xa,ya];

ra = ra.';
xflam = repelem(ra,1,2);


%scatter(ra(:,1),ra(:,2));
v0   = 1;
dcorr = 1/(4*pi)*log(dk^2*dlam^2/abs(4*pi^2- dk^2*dlam^2))+1i/4;
dcorr = dcorr*dk;

%%% now generate the plane wave...

vplane = exp(1i*dk*xa)*[0,1];


matfun = @(i,j)  afun_wrapper(i,j,ra,dk,v0,dcorr);
p = 8;
[pr,~,~,pin] = proxy_square_pts(4*p);
opdims = [2,2];
pxyfun = @(x,slf,nbr,l,ctr) proxyfun(xflam,slf,nbr,l,ctr,ra,dk,v0,opdims,pr,pin);

vplane2 = zeros(2*natoms,1);
vplane2(1:2:end) = vplane(1:natoms);
vplane2(2:2:end) = vplane(natoms+1:end);

occ = 100;
rank_or_tol = 1e-6;
opts = [];
opts.verb = 1;

tic, F = rskelf(matfun,xflam,occ,rank_or_tol,pxyfun,opts); toc
vsol2 = rskelf_sv(F,vplane2);
% Mnew = afun(1:natoms,1:natoms,ra,dk,v0,dcorr);
% vsol2 = Mnew\vplane2;

vsol_final = [vsol2(1:2:end) vsol2(2:2:end)];
[xss,yss] = meshgrid(xs,ys);
xss = xss(:);
yss = yss(:);
targs = [xss(:).';yss(:).'];
eps = rank_or_tol;
zk = complex(dk);
srcinfo = [];
srcinfo.nd = 2;
srcinfo.sources = ra;
srcinfo.charges = vsol_final.';
pg = 0 ;
pgt = 2;
[U] = hfmm2d(eps,zk,srcinfo,pg,targs,pgt);

vf1 = dk*squeeze(U.pottarg(1,:).') + 1j*squeeze((U.gradtarg(1,1,:) + U.gradtarg(2,2,:)));
vf2 = dk*squeeze(U.pottarg(2,:).') + 1j*squeeze((U.gradtarg(1,2,:) - U.gradtarg(2,1,:)));


vf1 = vf1*v0;
vf2 = vf2*v0;


[xss,yss] = meshgrid(xs,ys);
vplot = reshape(abs(vf1),size(xss));

h = pcolor(xss,yss,vplot);
set(h,'EdgeColor','none');    
hold on;
scatter(xa,ya,10,'filled');



vplane = exp(1i*dk*xss(:));
vf1_total = vf1;
vf2_total = vf2+vplane;

vplot = reshape(abs(vf1_total),size(xss));
figure(1);
clf
h = pcolor(xss,yss,vplot);
%caxis([0,4])
set(h,'EdgeColor','none');    
hold on;
scatter(xa,ya,5,'filled','white');
axis equal
colorbar()


vplot = reshape(abs(vf2_total),size(xss));
figure(2);
clf
h = pcolor(xss,yss,vplot);
%caxis([0,4])
set(h,'EdgeColor','none');    
hold on;
scatter(xa,ya,5,'filled','white');
axis equal
colorbar()

vplot = reshape(abs(vf1_total).^2 + abs(vf2_total).^2,size(xss));
figure(3);
clf();
h = pcolor(xss,yss,vplot);
caxis([0,4])
axis equal
xlim([xmin,xmax])
ylim([ymin,ymax])
set(h,'EdgeColor','none');    
hold on;
scatter(xa,ya,5,'filled','white');
colorbar()
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18);
saveas(gcf,join(['./results/reg_nx' int2str(nx0) '_ny' ...
     int2str(ny0) '_nover' int2str(nover) '_xdir_msign.pdf']));
fname_save = ['./data/everything_reg_nx' int2str(nx0) '_ny' int2str(ny0) ...
      '_nover' int2str(nover) '_xdir_msign.mat'];
save(fname_save);


