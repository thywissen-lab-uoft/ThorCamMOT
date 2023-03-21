if doBoxCount
    boxPopts = struct;
    boxPopts.FigLabel = FigLabel;
    boxPopts.xUnit=thor_unit;
    boxPopts.NumberExpFit = 0;        % Fit exponential decay to atom number
    boxPopts.NumberLorentzianFit=0;   % Fit atom number to lorentzian
    boxPopts.CenterSineFit = 0;       % Fit sine fit to cloud center
    boxPopts.CenterDecaySineFit = 0;  % Fit decaying sine to cloud center
    boxPopts.CenterParabolaFit = 0;
    boxPopts.CenterLinearFit = 0;     % Linear fit to cloud center
    boxPopts.NumberExpOffsetFit = 0; % Exp decay fit with nonzero offset    
       
    hF_number_box = showCounts(box_data,thor_xVar,boxPopts);  
    ylim([0 max(get(gca,'YLim'))]);    
    if doSave;saveFigure(hF_number_box,'box_number',saveOpts);end
    
    if ~isequal(thor_xVar,'ExecutionDate')
        hF_number_box_time = showCounts(...
            chDataXVar(box_data,'ExecutionDate'),'ExecutionDate',boxPopts);  
        hF_number_box_time.Position(2)=700;
        ylim([0 max(get(gca,'YLim'))]);    
        if doSave;saveFigure(hF_number_box_time,'box_number_time',saveOpts);end
    end
        
end
%%
if doGaussFit
    
    gaussPopts = struct;
    gaussPopts.FigLabel = FigLabel;
    gaussPopts.xUnit=thor_unit;
    gaussPopts.NumberExpFit = 0;        % Fit exponential decay to atom number
    gaussPopts.NumberLorentzianFit=0;   % Fit atom number to lorentzian
    gaussPopts.CenterSineFit = 0;       % Fit sine fit to cloud center
    gaussPopts.CenterDecaySineFit = 0;  % Fit decaying sine to cloud center
    gaussPopts.CenterParabolaFit = 0;
    gaussPopts.CenterLinearFit = 0;     % Linear fit to cloud center
    gaussPopts.NumberExpOffsetFit = 0; % Exp decay fit with nonzero offset    
       
    hF_number_gauss = showCounts(gauss_data,thor_xVar,gaussPopts);  
    ylim([0 max(get(gca,'YLim'))]);    
    if doSave;saveFigure(hF_number_gauss,'gauss_number',saveOpts);end
    
    if ~isequal(thor_xVar,'ExecutionDate')
        hF_number_gauss_time = showCounts(...
            chDataXVar(box_data,'ExecutionDate'),'ExecutionDate',gaussPopts);  
        hF_number_gauss_time.Position(2)=700;
        ylim([0 max(get(gca,'YLim'))]);    
        if doSave;saveFigure(hF_number_gauss_time,'gauss_number_time',saveOpts);end
    end
    
    [hF_center_gauss,hF_center_gaussb] = showCenter(gauss_data,thor_xVar,gaussPopts);  
    if doSave;saveFigure(hF_center_gauss,'gauss_center',saveOpts);end
    if doSave;saveFigure(hF_center_gaussb,'gauss_center2',saveOpts);end

    [hF_size_gauss] = showCenter(gauss_data,thor_xVar,gaussPopts);  
    if doSave;saveFigure(hF_size_gauss,'gauss_size',saveOpts);end
        
end