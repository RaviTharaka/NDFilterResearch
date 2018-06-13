function [] = visualize_3d(head,SIG,thresh,varargin)
    
    f = figure('name', head, 'position', [800 70 700 700]);
    w = (-1:2/size(SIG,1):1-2/size(SIG,1));   
    
    SIG_db = mag2db(abs(SIG));
    
    t = thresh;
    SIG_db_temp = SIG_db > t;
    
    [x,y,ct] = ind2sub(size(SIG_db_temp), find(SIG_db_temp));
    x = x/size(SIG,1)*2-1;
    y = y/size(SIG,2)*2-1;
    ct = ct/size(SIG,3)*2-1;
    
    plot3(x, y, ct, 'k.');

%     p = patch(isosurface(w,w,w,SIG_db,thresh));  
%     set(p,'FaceColor','blue','EdgeColor','none');
%     camlight; lighting phong;

    daspect([1 1 1]);
    grid on;
    axis ([-1,1,-1,1,-1,1]);
    xlabel('\omega_x (\times\pi), rad/sample');
    ylabel('\omega_y (\times\pi), rad/sample');
    zlabel('\omega_{ct} (\times\pi), rad/sample');
        
    text(1,1,1,strcat(num2str(t),'dB'));
    
    b = uicontrol('Parent',f,'Style','slider','Position',[300,10,200,25],...
              'value',t, 'min',-53, 'max',0 ...
              , 'callback', {@update,SIG});      
end


function update(es,ed,SIG,varargin)
    SIG_db = mag2db(abs(SIG));
    SIG_db_temp = SIG_db > es.Value;
    [x, y, ct] = ind2sub(size(SIG_db_temp), find(SIG_db_temp));
    x = x/size(SIG,1)*2-1;
    y = y/size(SIG,2)*2-1;
    ct = ct/size(SIG,3)*2-1;
    plot3(x, y, ct, 'k.');

    daspect([1 1 1]);
    grid on;
    axis ([-1,1,-1,1,-1,1]);
    xlabel('\omega_x (\times\pi), rad/sample');
    ylabel('\omega_y (\times\pi), rad/sample');
    zlabel('\omega_{ct} (\times\pi), rad/sample');
    
    text(1,1,1,strcat(num2str(es.Value),'dB'));
%     subplot(2,2,2);
%     pcolor(-1:2/size(s,1):1-1/size(s,1),-1:2/size(s,1):1-1/size(s,1),abs(s));
%     refline(0,round(es.Value)/size(s,1)*2-1);
%     axis square;
%     shading interp;
%     title 'Signal Fourier Domain';
%     xlabel '\omega_x (\times\pi)';
%     ylabel '\omega_t (\times\pi)';
%     
%     subplot(2,2,4);
%     plot(-1:2/size(s,1):1-1/size(s,1),mag2db(abs(s(round(es.Value),:))));
%     title 'Cross-section Fourier Domain';
%     xlabel '\omega_x (\times\pi)';
%     ylabel 'Magnitude (dB)';
%     axis([-1, 1, -70, 30]);
%     
%     subplot(2,2,3);
%     plot(-1:2/size(s,1):1-1/size(s,1),mag2db(abs(s(round(es.Value),:))));
%     title 'Passband Cross-section Fourier Domain';
%     xlabel '\omega_x (\times\pi)';
%     ylabel 'Magnitude (dB)';
%     if length(varargin)==2
%         axis([varargin{1} varargin{2} -5 5]);
%     else
%         axis([-1, 1, -5, 5]);
%     end
end

