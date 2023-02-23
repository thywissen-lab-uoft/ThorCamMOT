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
%     ylim([0 max(get(gca,'YLim'))]);    
    if doSave;saveFigure(hF_number_box,'box_number',saveOpts);end
    
    if ~isequal(thor_xVar,'ExecutionDate')
        hF_number_box_time = showCounts(...
            chDataXVar(box_data,'ExecutionDate'),'ExecutionDate',boxPopts);  
        hF_number_box_time.Position(2)=700;
%         ylim([0 max(get(gca,'YLim'))]);    
        if doSave;saveFigure(hF_number_box_time,'box_number_time',saveOpts);end
    end
        
end