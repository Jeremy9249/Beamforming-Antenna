%GedoraError1Interference
%Beamstearing using Gedora Version 1.00
%Accurary Comparison
%Last modified 4/7/19 -- Jeremy Perez

%Calculating Weights
theta=-90:90;

% Degree Representation Vector
size=fix(180/spacing)-1;
u = [1 0]';             % 1= Desired, 0=Interference

thetaI1V=zeros(1,size,1,1);
maxval=zeros(1,size,1,1);
index=zeros(1,size,1,1);
max2val=zeros(1,size,1,1);
index2=zeros(1,size,1,1);
null=zeros(1,size,1,1);
raderror=zeros(1,size,1,1);
raderrorg=zeros(1,size,1,1);

%Desired Interupts
thetaI1=-89;
%thetaI2=89;
for a = 1:size
% Create matrix of steering vectors 
A=[];
thetaS = [0 thetaI1];
for inde=1:length(u)
A=[A,ARV(thetaS(inde)+91,:).'];
end
% Determine Array Weights 
w = u.'*A'*inv(A*A'+1e-9*eye(4));  % add small noise to diagonal
                                   % of A*A' so inverse is non-singular                                 
                                   
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
hold on
end

plot(thetaI1V,raderrorg,'-','DisplayName','Godara','Color',colorstr{3})
xlabel('Angle (deg)')
ylabel('Error (dB)')
%title('Desired Null vs Response')
axis([-90 90 -30 0])
grid on
box on