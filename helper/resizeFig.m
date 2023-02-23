function resizeFig(hF,t,axs)

if nargin ~=3
    axs=[];    
end

function chSize(~,~)
    try
        t.Position(3)=t.Parent.Position(3);
        t.Position(4)=t.Extent(4);
        t.Position(1:2)=[5 t.Parent.Position(4)-t.Position(4)];
        drawnow;
        
       
        
        for kk=1:length(axs)
            axs(kk).Units='pixels';
            axs(kk).Position(2)=45;
            axs(kk).Position(4)=...
                (axs(kk).Parent.Position(4) - ...
                t.Position(4))-axs(kk).Position(2)-20;
            axs(kk).Units='normalized'; 
        end
    catch ME
        warning('resize issue');
    end 
end
chSize;

hF.SizeChangedFcn=@chSize;
end

