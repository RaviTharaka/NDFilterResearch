function [] = visualize(head,sig,SIG, varargin)
    SIG = abs(SIG);
    
    f = figure('name', head, 'position', [70 70 900 730]);
    subplot(2,2,1);
    surf(abs(sig));
    title 'Signal Time Domain';
    xlabel 'x';
    ylabel 't';
    zlabel 'Mag';
    axis tight;
    shading interp;
    
    subplot(2,2,2);
    pcolor(-1:2/size(SIG,2):1-1/size(SIG,2),-1:2/size(SIG,1):1-1/size(SIG,1),SIG);
    refline(0,-0.5);
    title 'Signal Fourier Domain';
    xlabel '\omega_x (\times\pi)';
    ylabel '\omega_t (\times\pi)';
    axis square;
    shading interp;
    
    subplot(2,2,3);
    plot(-1:2/size(SIG,2):1-1/size(SIG,2),mag2db(SIG(255,:)));
    title 'Passband Cross-section Fourier Domain';
    xlabel '\omega_x (\times\pi)';
    ylabel 'Magnitude (dB)';
    if length(varargin)==2
        axis([varargin{1} varargin{2} -5 5]);
    else
        axis([-1, 1, -5, 5]);
    end
    
    subplot(2,2,4);
    plot(-1:2/size(SIG,2):1-1/size(SIG,2),mag2db(SIG(255,:)));
    title 'Cross-section Fourier Domain';
    xlabel '\omega_x (\times\pi)';
    ylabel 'Magnitude (dB)';
    axis([-1, 1, -70, 30]);
    
    subplot(2,2,3);
        
    if length(varargin)==2
        b = uicontrol('Parent',f,'Style','slider','Position',[600,10,430,25],...
              'value',255, 'min',1, 'max',size(SIG,1), 'sliderstep',[1/size(SIG,1) 1/30]...
              , 'callback', {@update,SIG,varargin{1},varargin{2}});
    else
        b = uicontrol('Parent',f,'Style','slider','Position',[600,10,430,25],...
              'value',255, 'min',1, 'max',size(SIG,1), 'sliderstep',[1/size(SIG,1) 1/30]...
              , 'callback', {@update,SIG});
    end        
end

function update(es,ed,s,varargin)
    subplot(2,2,2);
    pcolor(-1:2/size(s,2):1-1/size(s,2),-1:2/size(s,1):1-1/size(s,1),s);
    refline(0,round(es.Value)/size(s,1)*2-1);
    axis square;
    shading interp;
    title 'Signal Fourier Domain';
    xlabel '\omega_x (\times\pi)';
    ylabel '\omega_t (\times\pi)';
    
    subplot(2,2,4);
    plot(-1:2/size(s,2):1-1/size(s,2),mag2db(s(round(es.Value),:)));
    title 'Cross-section Fourier Domain';
    xlabel '\omega_x (\times\pi)';
    ylabel 'Magnitude (dB)';
    axis([-1, 1, -70, 30]);
    
    subplot(2,2,3);
    plot(-1:2/size(s,2):1-1/size(s,2),mag2db(s(round(es.Value),:)));
    title 'Passband Cross-section Fourier Domain';
    xlabel '\omega_x (\times\pi)';
    ylabel 'Magnitude (dB)';
    if length(varargin)==2
        axis([varargin{1} varargin{2} -5 5]);
    else
        axis([-1, 1, -5, 5]);
    end
end

