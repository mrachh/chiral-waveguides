ys = -20:0.02:20;
dh = ys(2)-ys(1);

sig= 0.5;

vg = exp(-ys.^2/(2*sig^2));
de = 0.1;


v1 = vg;
v2 = 0*vg;

v10 = v1;
v20 = v2;

dt = dh/256;
nt = 5/dt;


for ii=1:nt
    [rhs1,rhs2] = get_rhs_11(v1,v2,de,dh);
    v1t = v1 + dt*rhs1;
    v2t = v2 + dt*rhs2;
    for ijk=1:4
    [rhs1t,rhs2t] = get_rhs_11(v1t,v2t,de,dh);
    v1t = v1 + dt*(rhs1+rhs1t)/2;
    v2t = v2 + dt*(rhs2+rhs2t)/2;
    end
    v1 = v1t;
    v2 = v2t;
end

ttot = nt*dt;

[jf11,jf21,jf12,jf22] = gf_pde_smth(ys,ttot,de);
jf11 = circshift(jf11,find(ys==0));
jf21 = circshift(jf21,find(ys==0));
jf12 = circshift(jf12,find(ys==0));
jf22 = circshift(jf22,find(ys==0));

vfin1 = ifft(fft(jf11).*fft(v10)+fft(jf12).*fft(v20))*dh;
vfin2 = ifft(fft(jf21).*fft(v10)+fft(jf22).*fft(v20))*dh;

figure(1)
clf
plot(ys,real(vfin1-v1),ys,imag(vfin1-v1),ys,real(vfin2-v2),ys,imag(vfin2-v2)); 

vdiff1 = v1-vfin1;
vdiff2 = v2-vfin2;

%%%%% now post-process
vdiffl1 = circshift(vdiff1,find(ys==ttot)-find(ys==0));
vdiffr1 = circshift(vdiff1,find(ys==-ttot)-find(ys==0));

vdiffl2 = circshift(vdiff2,find(ys==ttot)-find(ys==0));
vdiffr2 = circshift(vdiff2,find(ys==-ttot)-find(ys==0));

vdiffl1_wg = vdiffl1 - v10/2*exp(-1j*de*ttot) + v20/2*exp(-1j*de*ttot);
vdiffl2_wg = vdiffl2 + v10/2*exp(-1j*de*ttot) - v20/2*exp(-1j*de*ttot);

vdiffr1_wg = vdiffr1 - v10/2*exp(-1j*de*ttot) - v20/2*exp(-1j*de*ttot);
vdiffr2_wg = vdiffr2 - v10/2*exp(-1j*de*ttot) - v20/2*exp(-1j*de*ttot);


figure(2)
clf
plot(ys,real(vdiffl1),'k.'); hold on; plot(ys,imag(vdiffl1),'r.');
plot(ys,real(vdiffl1_wg),'k--'); plot(ys,imag(vdiffl1_wg),'r--'); 


figure(3)
clf
plot(ys,real(vdiffl2),'k.'); hold on; plot(ys,imag(vdiffl2),'r.');
plot(ys,real(vdiffl2_wg),'k--'); plot(ys,imag(vdiffl2_wg),'r--'); 


figure(4)
clf
plot(ys,real(vdiffr1),'k.'); hold on; plot(ys,imag(vdiffr1),'r.');
plot(ys,real(vdiffr1_wg),'k--'); plot(ys,imag(vdiffr1_wg),'r--'); 


figure(5)
clf
plot(ys,real(vdiffr2),'k.'); hold on; plot(ys,imag(vdiffr2),'r.');
plot(ys,real(vdiffr2_wg),'k--'); plot(ys,imag(vdiffr2_wg),'r--'); 

emaxl1 = max(abs(vdiffl1(abs(ys)<=2)));
emaxl1_wg =  max(abs(vdiffl1_wg(abs(ys)<=2)));

emaxl2 = max(abs(vdiffl2(abs(ys)<=2)));
emaxl2_wg =  max(abs(vdiffl2_wg(abs(ys)<=2)));

emaxr1 = max(abs(vdiffr1(abs(ys)<=2)));
emaxr1_wg =  max(abs(vdiffr1_wg(abs(ys)<=2)));

emaxr2 = max(abs(vdiffr2(abs(ys)<=2)));
emaxr2_wg =  max(abs(vdiffr2_wg(abs(ys)<=2)));


errs = [emaxl1, emaxl2; emaxr1, emaxr2];
errs_wg = [emaxl1_wg, emaxl2_wg; emaxr1_wg, emaxr2_wg];
