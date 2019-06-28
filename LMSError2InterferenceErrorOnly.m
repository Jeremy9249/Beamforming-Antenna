%BeamSteering2 using LMS
%Comparison of 2 Interrupts
%Last Modified 4/7/19 Jeremy Perez

%Calculating Weights
%N = input(' How many element do you want in uniform linear array? '); % number of elements in array
N= 4;
thetaS = 0;
size=fix(180/spacing);

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
%----- Desired Signal & Interferer -----%
T=1E-3;
t=(1:100)*T/100;
it=1:100;
S=cos(2*pi*t/T);
I = randn(1,100); 

%----- Create Array Factors for each user's signal for linear array -----%

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

e = conj(S(n)) - y; esARVe(n) = abs(e)^2;
% w = w +mu*e*conj(x);
w=w+mu(n)*conj(e)*x;
wi(:,n)=w;
yy(n)=y;

mu(n);
x1 = S(n)*vS ;% + I(n)*vI;
%y = w*x.';
y1=w'*x;

e1 = conj(S(n)) - y1;
esARVe(n) = abs(e1)^2;
%w = w +mu*e*conj(x);
w1=w+mu(n)*conj(e1)*x1;
wi(:,n)=w1;
yy1(n)=y1;

end
w = (w./w(1));% normalize results to first weight
w1 = (w1./w1(1));
%----- Plot Results -----%

theta = -pi/2:.02:pi/2;
AF = zeros(1,length(theta));
% Determine the array factor for linear array

AF = w(4)'.*ant1_cal+w(3)'.*ant2_cal+w(2)'.*ant3_cal+w(1)'.*ant4_cal;
theta=-90:90;

[maxval(a),index(a)] = max(abs(AF)/max(abs(AF)));
[pks,locs] = findpeaks(abs(AF)/max(abs(AF)));
pklocs=[pks,locs];
spklocs=sortrows(pklocs,'descend'); 
index(a)=index(a)-91;

responsev=20*log10(abs(AF)/max(abs(AF)));
thetaI1V(a)=thetaI1;
raderror(a)=responsev(thetaI1+91);
if raderror(a)<-30
    raderrorg(a)=-30;
else
    raderrorg(a)=raderror(a);
end

thetaI2V(a)=thetaI2;
rad2error(a)=responsev(thetaI2+91);
if rad2error(a)<-30
    rad2errorg(a)=-30;
else
    rad2errorg(a)=rad2error(a);
end

thetaI1=thetaI1+spacing;
thetaI2=thetaI2-spacing;
hold on
end
plot(thetaI1V,raderrorg,'-','DisplayName','LMS','Color',colorstr{1})
xlabel('Angle (deg)')
ylabel('Error (dB)')
% title('Desired Null 1 vs Response')
axis([-90 90 -30 0])
grid on
hold on