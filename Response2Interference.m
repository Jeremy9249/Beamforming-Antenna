%BeamSteering2 using LMS
%Comparison of 1 Interrupt
%Last Modified 4/7/19 Jeremy Perez

%Calculating Weights
%N = input(' How many element do you want in uniform linear array? '); % number of elements in array
N= 4;
u=[1 0 0]';
thetaS = 0;

%Desired Interupts
thetaI1 = -60;
thetaI2 =  50;
fig=figure;

for a = 1:h%:size+1
    % Create matrix of steering vectors
    
    %----- Desired Signal & Interferer -----%
    T=1E-3;
    t=(1:100)*T/100;
    it=1:100;
    S=cos(2*pi*t/T);
    I = randn(1,100);
    
    %----- Create Array Factors for each user's signal for linear array -----%
    vS = []; vI = [];vI1=[];
    vS=ARV(thetaS+91,:).';
    vI1=ARV(thetaI1+91,:).';
    vI2=ARV(thetaI2+91,:).';
    vI = vI1+vI2;
    
    %----- Solve for Weights using LMS -----%
    w = zeros(N,1);
    snr = 20e-6; % signal to noise ratio
    X=(vS+vI);   % Vector of desired and undesired angles
    %X=vS;   %no interupts
    Rx=X*X';     % Matrix of the vetor times its tranpose?
    mu=1/(real(trace(Rx)));
    
    wi=zeros(N,max(it));
    oldmu = mu;
    for n = 1:100
        mu(n) = oldmu/(1-(oldmu^(n+1)));
        oldmu = mu(n);
    end
    
    for n = 1:length(S)
        x = S(n)*vS + I(n)*vI;
        %y = w*x.';
        y=w'*x;
        
        et = conj(S(n)) - y;
        % w = w +mu*e*conj(x);
        w=w+mu(n)*conj(et)*x;
        wi(:,n)=w;
        mu(n);
        x1 = S(n)*vS ;% + I(n)*vI;
        %y = w*x.';
        y1=w'*x;
        
        e1 = conj(S(n)) - y1;
        esARVe(n) = abs(e1)^2;
        %w = w +mu*e*conj(x);
        w1=w+mu(n)*conj(e1)*x1;
        wi(:,n)=w1;
    end
    w = (w./w(1));% normalize results to first weight
    w1 = (w1./w1(1));
    %----- Plot Results -----%
    theta = -pi/2:.02:pi/2;
    AF = zeros(1,length(theta));
    % Determine the array factor for linear array
    AF = w(4)'.*ant1_cal+w(3)'.*ant2_cal+w(2)'.*ant3_cal+w(1)'.*ant4_cal;
    theta=-90:90;
    legendstr=sprintf('LMS', thetaI1);
    plot(theta,20*log10(abs(AF)/max(abs(AF))),'Color',colorstr{1}...
        ,'DisplayName',legendstr,'LineWidth',2)
    %legend('Location','northoutside')
    set(gca,'FontSize',12)
    set(gca,'xtick',[-90 -60 -30 0 30 60 90])
    xlabel('Angle (deg)')
    ylabel('Normalized Radiation Pattern (dB)')
    axis([-90 90 -30 0])
    %title('Array Factor Response 1 Null')
    %set(gca,'xtick',[-90 -60 -30 0 30 60 90])
    grid on
    hold on
    
    A=[];
    thetaSV = [0 thetaI1 thetaI2];
    for inde=1:length(u)
        A=[A,ARV(thetaSV(inde)+91,:).'];
    end
    % Determine Array Weights
    w = u.'*A'*inv(A*A'+1e-9*eye(4));  % add small noise to diagonal
    % of A*A' so inverse is non-singular
    
    % Determine Array Factor
    y = w(1)'.*ant1_cal + w(2)'.*ant2_cal + w(3)'.*ant3_cal + w(4)'.*ant4_cal;
    
    % Plot Results
    legendstr=sprintf('Godara', thetaI1);
    plot(theta,20*log10(abs(y)/max(abs(y))),'DisplayName',...
        legendstr,'Color',colorstr{3},'LineWidth',1.5)
    %title(letterstr{h},'FontWeight','normal');
    %legend('Location','northoutside')
    set(gca,'FontSize',12)
    set(gca,'xtick',[-90 -60 -30 0 30 60 90])
    axis([-90 90 -30 0])
    %title('\bfArray Factor Pattern for Fixed Weight Beamformer')
    xlabel('Angle (deg)');
    ylabel('Normalized Pattern (dB)')
    grid on
    hold off
    legend('Location','northoutside','Orientation','horizontal','NumColumnsMode','manual','NumColumns',2)
    set(gcf, 'color', 'white');

    drawnow
    frame=getframe(fig);
    im{a}=frame2im(frame);
    
    thetaI1=thetaI1+10;
    thetaI2=thetaI2-10;
end

close;
figure;
for idx=1:h
        subplot(5,5,idx)
        imshow(im{idx});
end

filename = 'TwoIntScan.gif'; % Specify the output file name
for idx = 1:h
    [AAA,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(AAA,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(AAA,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end
