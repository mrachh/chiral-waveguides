clear

A = load('hyperbolic/data/atoms.mat');
xs = 0:0.25/4:15;
ys = -25:0.25/4:25;
natoms = A.natoms;
xa = A.xatoms;
ya = A.yatoms;
[xasort,isort] = sort(xa);
ra = [xa,ya];

ra = ra.';
xflam = repelem(ra,1,2);


%scatter(ra(:,1),ra(:,2));

dlam = A.s0;
dk   = A.de;
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

occ = 200;
rank_or_tol = 1e-6;
opts = [];
opts.verb = 1;

tic, F = rskelf(matfun,xflam,occ,rank_or_tol,[],opts); toc
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
vf2_total = vf2 + vplane;

vplot = reshape(abs(vf1_total),size(xss));
figure;
h = pcolor(xss,yss,vplot);
caxis([0,4])
set(h,'EdgeColor','none');    
hold on;
scatter(xa,ya,40,'filled','white');
axis equal
colorbar()


vplot = reshape(abs(vf2_total),size(xss));
figure;
h = pcolor(xss,yss,vplot);
caxis([0,4])
set(h,'EdgeColor','none');    
hold on;
scatter(xa,ya,40,'filled','white');
axis equal
colorbar()

vplot = reshape(abs(vf1_total).^2 + abs(vf2_total).^2,size(xss));
figure;
h = pcolor(xss,yss,vplot);
caxis([0,4])
axis equal
xlim([min(xs(:)),max(xs(:))]);
ylim([min(ys(:)),max(ys(:))]);
set(h,'EdgeColor','none');    
hold on;
scatter(xa,ya,40,'filled','white');
colorbar()
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',18);

saveas(gcf,'./results/ref_elliptic_xdir_msign.pdf');
     



