%Beamstearing using Both Methods 
%Accurary Comparison
%Last modified 4/10/19 -- Jeremy Perez

%Calculating Weights
theta=-90:90;
N=4;
thetaS = [-80];   % Degree Representation Vector
%size=fix(180/spacing);
u = [1]';             % 1= Desired, 0=Interference
fig=figure;
for a = 1:h%:size+1
% Create matrix of steering vectors 
A=[];
for inde=1:length(u)
A=[A,ARV(thetaS(inde)+91,:).'];
end
% Determine Array Weights 
w = u.'*A'*inv(A*A'+1e-9*eye(4));  % add small noise to diagonal
                                   % of A*A' so inverse is non-singular                                  
disp(['  The array weights for Godara ' ,num2str(thetaS), ' beamformer are:'])
disp(' ')
for m = 1:length(w)
    disp(['   w',num2str(m),' = Mag:',num2str(abs(w(m))),' Phase:',num2str(rad2deg(angle(w(m))))])
end

% Determine Array Factor
y = w(1)'.*ant1_cal + w(2)'.*ant2_cal + w(3)'.*ant3_cal + w(4)'.*ant4_cal;

% Plot Results 
legendstr=sprintf('Godora', thetaS);
%subplot(2,2,a+p)
plot(theta,(20*log10(abs(y)/max(abs(y)))),'Color',colorstr{3},'DisplayName',legendstr)
xlabel({'Angle (deg)'})
axis([-90 90 -30 0])
xticks([-90 -60 -30 0 30 60 90])
xlabel('Angle (deg)');
ylabel('Normalized Pattern (dB)')
grid on
%legend('Location','northoutside')
set(gca,'FontSize',12)
%----- Desired Signal & Interferer -----%
T=1E-3;
t=(1:100)*T/100;
it=1:100;
S=cos(2*pi*t/T);
I = randn(1,100); 

%----- Create Array Factors for each user's signal for linear array -----%
% vI = [];
vS=ARV(thetaS+91,:).';

%----- Solve for Weights using LMS -----%
w = zeros(N,1); 
snr = 20e-6; % signal to noise ratio
X=vS;   %no interupts
Rx=X*X';     % Matrix of the vetor times its tranpose?
mu=1/(real(trace(Rx)));
oldmu = mu;
for n = 1:100
mu(n) = oldmu/(1-(oldmu^(n+1)));
oldmu = mu(n);
end

for n = 1:length(S)
x = S(n)*vS ;%+ I(n)*vI;
y=w'*x;
e = conj(S(n)) - y;
w=w+mu(n)*conj(e)*x;

end
w = (w./w(1));% normalize results to first weight
%disp(['  The array weights for LMS ' ,thetaS,' beamformer are:'])
%disp(' ')
% for m = 1:length(w)
%     disp(['   w',num2str(m),' = ',num2str(conj(w(m)),3)])
% end

%----- Plot Results -----%
% Determine the array factor for linear array
AF = w(4)'.*ant1_cal+w(3)'.*ant2_cal+w(2)'.*ant3_cal+w(1)'.*ant4_cal;
hold on
theta=-90:90;
legendstr=sprintf('LMS', thetaS);
plot(theta,(20*log10(abs(AF)/max(abs(AF)))),'Color',colorstr{1},...
     'DisplayName',legendstr,'LineWidth',1.5)
xlabel({'Angle (deg)'})
ylabel('Normalized Pattern (dB)')
axis([-90 90 -30 0])
xticks([-90 -60 -30 0 30 60 90])
yticks([-30 -20 -10 0])
grid on
set(gca,'FontSize',12)
hold off
legend('Location','northoutside','Orientation','horizontal','NumColumnsMode','manual','NumColumns',2)
set(gcf, 'color', 'white');
drawnow
frame=getframe(fig);
im{a}=frame2im(frame);

thetaS=thetaS+10;

end
close;
figure;
for idx=1:h
        subplot(5,5,idx)
        imshow(im{idx});
end

filename = 'BeamScan.gif'; % Specify the output file name
for idx = 1:h
    [AAA,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(AAA,map,filename,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(AAA,map,filename,'gif','WriteMode','append','DelayTime',1);
    end
end