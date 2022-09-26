ys = -20:0.005:20;
dh = ys(2)-ys(1);

sig= 0.5;

vg = exp(-ys.^2/(2*sig^2));
de = 0.1;


x0 = 5;
y0 = 0;
s0 = 0.1/2;
v0 = 1/2/pi/s0/s0;

th = 0;

v1 = 0*vg;
v2 = 0*vg;

v10 = v1;
v20 = v2;

dt = dh/2;
nt = 15/dt;

vals1=[];
vals2=[];
ts   =[];

t = 0;

for ii=1:nt
    [rhs1,rhs2] = get_rhs_11(v1,v2,de,dh);
    [v1src,v2src] = get_pde_src(t,ys,x0,y0,s0,v0,v1,v2,de,th);
    v1t = v1 + dt*(rhs1+v1src);
    v2t = v2 + dt*(rhs2+v2src);
    for ijk=1:2
    tt = t + dt;    
    [rhs1t,rhs2t] = get_rhs_11(v1t,v2t,de,dh);
    [v1srct,v2srct] = get_pde_src(tt,ys,x0,y0,s0,v0,v1t,v2t,de,th);
    v1t = v1 + dt*(rhs1+rhs1t+v1src+v1srct)/2;
    v2t = v2 + dt*(rhs2+rhs2t+v2src+v2srct)/2;
    end
    v1 = v1t;
    v2 = v2t;
    t = t + dt;
    if (mod(ii,round(nt/100)) == 0)
        vals1 = [vals1;v1];
        vals2 = [vals2;v2];
        ts    = [ts,t];
    end
end

ttot = nt*dt;

[vr1,vr2] = get_pde_src(x0,y0,x0,y0,s0,v0,0,0,de,th);
% Extract just the plane wave part of the data
vi1 = vr1/-1j/v0;
vi2 = vr2/-1j/v0;
delt = 0.001;
Gfun =  get_gfunself(s0,delt,de);
Gfun = Gfun/2/pi/2/pi; % to account for the leading 1/(2pi)^2 without slowing down adaptive integration
% M = eye(2) - Gfun;
% M2 = eye(2);
% M2(1,1) = M2(1,1) - Gfun(2,2);
% M2(2,2) = M2(2,2) - Gfun(1,1);
% 
% vr1 = vr1*(Gfun(1,1));
% vr2 = vr2*(Gfun(2,2));
% 
% vv = [vr1;vr2];
% sval = M\vv;

sval = complex(zeros(2,1));
sval(1) = -Gfun(1,1)*vi1/(1+Gfun(1,1)); %Implementing formula 8
sval(2) = -Gfun(2,2)*vi2/(1+Gfun(2,2));
sval(1) = sval(1) + vi1;
sval(2) = sval(2) + vi2;






% figure; h=pcolor(abs(vals1));set(h,'EdgeColor','none'); caxis([0,0.2]);
% figure; h=pcolor(abs(vals2));set(h,'EdgeColor','none'); caxis([0,0.2]);

%%%%%% generate the incoming field

    [XX,YY] = meshgrid(ts,ys);
    XX = XX.';
    YY = YY.';
    z = cos(th);
    k = sin(th);
    p = z-de;
    
    amp1 = sin(th);
    amp2 = 1-cos(th);
    ampn = sqrt(amp1^2+amp2^2);
    amp1 = amp1/ampn;
    amp2 = amp2/ampn;
    
    if (th == 0) 
        amp1 = 1;
        amp2 = 0;
    end
    
    v10 = amp1*exp(1i*p*XX+1i*k*YY);
    v20 = amp2*exp(1i*p*XX+1i*k*YY);
    
    
    XXscat = XX - x0;
    YYscat = YY - y0;
    rr = XXscat.^2 - YYscat.^2;
    
    XXscat(rr<=0) = 0;
    YYscat(rr<=0) = 0;
    rr(rr<=0) = 0;
    uu = sqrt(rr);
    
    
    J11=exp(-1j*de*XXscat)/2.*( ...
      1j*besselj(1,uu)./uu.*XXscat + besselj(0,uu));
    J22=exp(-1j*de*XXscat)/2.*( ...
      1j*besselj(1,uu)./uu.*XXscat - besselj(0,uu));
    J12=exp(-1j*de*XXscat)/2.*1j.*besselj(1,uu)./uu.*YYscat;
    J21=exp(-1j*de*XXscat)/2.*1j.*besselj(1,uu)./uu.*YYscat;
    
    J11(rr<=0) = 0;
    J12(rr<=0) = 0;
    J21(rr<=0) = 0;
    J22(rr<=0) = 0;
    
    
    vscat1 = sval(1)*J11 + sval(2)*J12;
    vscat2 = sval(1)*J21 + sval(2)*J22;
    vscat1 = vscat1;
    vscat2 = vscat2;
    vscat1(XXscat<=0) = 0;
    vscat2(XXscat<=0) = 0;
% figure; h=pcolor(abs(vscat1));set(h,'EdgeColor','none'); caxis([0,0.2]);
% figure; h=pcolor(abs(vscat2));set(h,'EdgeColor','none'); caxis([0,0.2]);
%      
% figure; h=pcolor(real(vscat1./vals1));set(h,'EdgeColor','none'); caxis([0,1]);
% figure; h=pcolor(real(vscat2./vals2));set(h,'EdgeColor','none'); caxis([0,1]);
%  
 
figure; h=pcolor(log10(abs(vscat1-vals1)));set(h,'EdgeColor','none'); caxis([-5,0]);
figure; h=pcolor(log10(abs(vscat2-vals2)));set(h,'EdgeColor','none'); caxis([-5,0]);


figure; plot(ys,real(vscat1(end,:)),'k.',ys,real(vals1(end,:)),'r.');
figure; plot(ys,imag(vscat1(end,:)),'k.',ys,imag(vals1(end,:)),'r.');
    
    %%%%% compute the total field
    
    vtot1 = vals1 + v10;
    vtot2 = vals2 + v20;
    
    vtot_g1 = vscat1 + v10;
    vtot_g2 = vscat2 + v20;
    
% figure; h=pcolor(abs(vtot1));set(h,'EdgeColor','none');
% figure; h=pcolor(abs(vtot2));set(h,'EdgeColor','none');
    
    ccur =abs(vtot1).^2+abs(vtot2).^2;
    
% figure; h=pcolor(abs(ccur));set(h,'EdgeColor','none');
%     
% figure; plot(sum(ccur.'));
% 

    
    ccur2 =abs(vtot_g1).^2+abs(vtot_g2).^2;
    

% figure; h=pcolor(abs(ccur2));set(h,'EdgeColor','none');
%     
%  figure; plot(sum(ccur2.'));
