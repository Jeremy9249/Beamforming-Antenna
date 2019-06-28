%Beamstearing using Gedora Version 1.00
%Accurary Comparison
%Last modified 4/7/19 -- Jeremy Perez

%Calculating Weights
theta=-90:90;

thetaS = [-89];   % Degree Representation Vector
size=fix(180/spacing)-1;
u = [1]';             % 1= Desired, 0=Interference

thetaSV=zeros(1,size,1,1);
maxval=zeros(1,11,1,1);
index=zeros(1,size,1,1);
max2val=zeros(1,size,1,1);
index2=zeros(1,size,1,1);
lobe=zeros(1,size,1,1);
perror=zeros(1,size,1,1);

for a = 1:size
% Create matrix of steering vectors 
A=[];
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
dpklocs=sortrows(pklocs,'descend'); %not right
index(a)=index(a)-91;
index2(a)=dpklocs(2,2)-91;
lobe(a)=dpklocs(2,1);
thetaSV(a)=thetaS;
perror(a)=abs(index(a)-thetaSV(a));
thetaS=thetaS+spacing;
hold on
end
plot(thetaSV,smoothdata(perror),'DisplayName','Godara','Color',colorstr{3},'Linewidth',1.5)
box on
xlabel('Angle (deg)')
ylabel('Error (deg)')
%title('Desired Angle vs Accuracy')
grid on
axis([-90 90 0 90])