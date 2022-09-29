iflam = 0;
if(exist('./FLAM/startup.m','file'))
    run './FLAM/startup.m'
    iflam = 1;
else
    str = which('rskelf');
    if(~isempty(str))
        iflam = 1;
    end
end

addpath(genpath('./src'));
    

