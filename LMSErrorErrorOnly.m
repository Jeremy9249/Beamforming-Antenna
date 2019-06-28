%BeamSteering using LMS   Version 1.0.0
%Accuracy Comparison
%Last Modified 3/27/19 Jeremy Perez
%Error Only

%Calculating Weights
%N = input(' How many element do you want in uniform linear array? '); % number of elements in array
N= 4;
thetaS  = -89;
size=fix(180/spacing)-1;

thetaSV=zeros(1,size,1,1);
maxval=zeros(1,size,1,1);
index=zeros(1,size,1,1);
max2val=zeros(1,size,1,1);
index2=zeros(1,size,1,1);
lobe=zeros(1,size,1,1);
perror=zeros(1,size,1,1);
for a = 1:size
%----- Desired Signal & Interferer -----%
T=1E-3;
t=(1:100)*T/100;
it=1:100;
S=cos(2*pi*t/T);
I = randn(1,100); 

%----- Create Array Factors for each user's signal for linear array -----%
vS = []; 
% vI = [];
vS=ARV(thetaS+91,:).';

%----- Solve for Weights using LMS -----%
w = zeros(N,1); 
snr = 20e-6; % signal to noise ratio
X=vS;   %no interupts
Rx=X*X';     % Matrix of the vetor times its tranpose?
mu=1/(real(trace(Rx)));

wi=zeros(N,max(it));
oldmu = mu;
for n = 1:100
mu(n) = oldmu/(1-(oldmu^(n+1)));

oldmu = mu(n);

end

for n = 1:length(S)
x = S(n)*vS ;%+ I(n)*vI;
y=w'*x;

e = conj(S(n)) - y;
esARVe(n) = abs(e)^2;
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

end
w = (w./w(1));% normalize results to first weight
w1 = (w1./w1(1));
%----- Plot Results -----%

% Determine the array factor for linear array
AF = w(4)'.*ant1_cal+w(3)'.*ant2_cal+w(2)'.*ant3_cal+w(1)'.*ant4_cal;

theta=-90:90;

[maxval(a),index(a)] = max(abs(AF)/max(abs(AF)));
[pks,locs] = findpeaks(abs(AF)/max(abs(AF)));
pklocs=[pks,locs];
dpklocs=sortrows(pklocs,'descend'); 
index(a)=index(a)-91;
lobe(a)=dpklocs(2,1);
thetaSV(a)=thetaS;
perror(a)=abs(index(a)-thetaSV(a));
thetaS=thetaS+spacing;
hold on
end

plot(thetaSV,smoothdata(perror),'DisplayName','LMS','Color',colorstr{1},'Linewidth',2)
xlabel('Angle (deg)')
ylabel('Error (deg)')
grid on
axis([-90 90 0 40])