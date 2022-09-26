function [v1t,v2t] = get_pde_src(x,ys,x0,y0,s0,v0,v1,v2,de,th)

%%%     plane wave properties
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
    
    v10 = amp1*exp(1i*p*x+1i*k*ys)+v1;
    v20 = amp2*exp(1i*p*x+1i*k*ys)+v2;
    
    xds = repmat((x-x0).^2,1,numel(ys));
    yds = (repmat(ys,numel(y0),1)-repmat(y0,1,numel(ys))).^2;
    vep = 1i*v0*exp(-(xds+yds)/(2*s0^2));
    v1s = sum(vep,1);
    v2s = sum(vep,1);
    %v1s = 1i*v0*exp(-((x-x0)^2+(ys-y0).^2)/(2*s0^2));
    %v2s = 1i*v0*exp(-((x-x0)^2+(ys-y0).^2)/(2*s0^2));

    v1t = v10.*v1s;
    v2t = v20.*v2s;

end

