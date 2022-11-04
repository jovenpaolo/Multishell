clear all;

f = waitbar(0,'1','Name','Wait laangs...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);
c=0;

% Update radFe
radFe=linspace( 1, 25, 26 );
lam=linspace( 300, 900, 601 )*1e-9;
radAu=linspace( 0, 25, 21 );

% Allocate matrices for speed
sct=zeros(length(radFe),length(radAu),length(lam));
abt=zeros(length(radFe),length(radAu),length(lam));
ext=zeros(length(radFe),length(radAu),length(lam));
Emax=zeros((length(radAu)*length(radFe)*length(lam)),4);

for i=1:length(radFe)
    for j=1:length(radAu)
        for k=1:length(lam)
        % Main Function
            [Emaxtable] = PeakE( radFe(i), radAu(j), lam(k));
        % Update Progress bar 
            c=c+1;
            progress=c/(length(radAu)*length(radFe)*length(lam));
        if getappdata(f,'canceling')
            break
        end
        % Know which part of the loop
        [i,j,k]

        Emax(c,1)=radFe(i);
        Emax(c,2)=radAu(j);
        Emax(c,3)=lam(k);
        Emax(c,4)=Emaxtable;

        waitbar(progress,f,sprintf('%12.2f',progress*100))
        end
    end
end

delete(f)
