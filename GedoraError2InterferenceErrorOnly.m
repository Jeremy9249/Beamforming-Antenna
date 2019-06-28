%BeamSteering2 using Gedora
%Comparison of 2 Interrupts
%Last Modified 4/7/19 Jeremy Perez

%Calculating Weights
%N = input(' How many element do you want in uniform linear array? '); % number of elements in array
N= 4;
thetaS = 0;
size=fix(180/spacing);
u = [1 0 0]';             % 1= Desired, 0=Interference

thetaI1V=zeros(1,size,1,1);
thetaI2V=zeros(1,size,1,1);
maxval=zeros(1,size,1,1);
index=zeros(1,size,1,1);
max2val=zeros(1,size,1,1);
index2=zeros(1,size,1,1);
null=zeros(1,size,1,1);
raderror=zeros(1,size,1,1);
raderrorg=zeros(1,size,1,1);
rad2error=zeros(1,size,1,1);
rad2errorg=zeros(1,size,1,1);
%Desired Interupts
thetaI1 = -89;
thetaI2 = 89;
for a = 1:size
% Create matrix of steering vectors 
A=[];
thetaS = [0 thetaI1 thetaI2];
for inde=1:length(u)
A=[A,ARV(thetaS(inde)+91,:).'];
end
% Determine Array Weights 
w = u.'*A'*inv(A*A'+1e-9*eye(4));  % add small noise to diagonal
                                   % of A*A' so inverse is non-singular                                 
                                   
% Display normalized weights 
% disp('  The array weights for the Maximum SIR beamformer are:')
% disp(' ')
% for m = 1:length(w)
%     disp(['   w',num2str(m),' = ',num2str(conj(w(m)))])
% end
% We display conj(w) because the weights solved for above are actually
% equal to w' where (') is the Hermitian transpose operator  

% Determine Array Factor
y = w(1)'.*ant1_cal + w(2)'.*ant2_cal + w(3)'.*ant3_cal + w(4)'.*ant4_cal;

[maxval(a),index(a)] = max(abs(y)/max(abs(y)));
[pks,locs] = findpeaks(abs(y)/max(abs(y)));
pklocs=[pks,locs];
spklocs=sortrows(pklocs,'descend'); 
index(a)=index(a)-91;

responsev=20*log10(abs(y)/max(abs(y)));
thetaI1V(a)=thetaI1;
raderror(a)=responsev(thetaI1+91);
if raderror(a)<-30
    raderrorg(a)=-30;
else
    raderrorg(a)=raderror(a);
end
thetaI1=thetaI1+spacing;
thetaI2=thetaI2-spacing;
hold on
end
plot(thetaI1V,raderrorg,'-','DisplayName','Godara','Color',colorstr{3})
xlabel('Angle (deg)')
ylabel('Error (dB)')
axis([-90 90 -30 0])
grid on
box on
% figure
% plot(thetaI2V,rad2errorg,'x')
% xlabel('AOI (deg)')
% ylabel('Response at AOI (dB)')
% title('Desired Null 2 vs Response')
% axis([-90 90 -30 0])
% grid on