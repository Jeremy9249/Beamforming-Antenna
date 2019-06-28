%Beam Steering Master File
clear,clc,close all;
%Read in IQ data from each of the channels
load('C:\Users\JP\Documents\MATLAB\Calibration\Captures1\ant1capture.mat')
load('C:\Users\JP\Documents\MATLAB\Calibration\Captures1\ant2capture.mat')
load('C:\Users\JP\Documents\MATLAB\Calibration\Captures1\ant3capture.mat')
load('C:\Users\JP\Documents\MATLAB\Calibration\Captures1\ant4capture.mat')
ant1=smoothdata(ant1);
ant2=smoothdata(ant2);
ant3=smoothdata(ant3);
ant4=smoothdata(ant4);

% Load variables
freq=2e9;
nsamps=length(ant1(1,:));
broadsideInd=91;
ARV=[];
% Find fourier transforms at broadside should be 0, it wont be
% ant1_FFT is the fourier vector of the antenna
% ant1 is the max amplitude of the vecotor
% ant1b is when it occurs
% ant_Phase returns the angle in degrees
ant1_FFT = fft(ant1(broadsideInd,:));
[ant1a,ant1b] = max(abs(ant1_FFT));
ant1_Phase = angle(ant1_FFT(ant1b))*180/pi;

ant2_FFT = fft(ant2(broadsideInd,:));
[ant2a,ant2b] = max(abs(ant2_FFT));
ant2_Phase = angle(ant2_FFT(ant2b))*180/pi;

ant3_FFT = fft(ant3(broadsideInd,:));
[ant3a,ant3b] = max(abs(ant3_FFT));
ant3_Phase = angle(ant3_FFT(ant3b))*180/pi;

ant4_FFT = fft(ant4(broadsideInd,:));
[ant4a,ant4b] = max(abs(ant4_FFT));
ant4_Phase = angle(ant4_FFT(ant4b))*180/pi;

% Calibrate Amplitudes at Broadside
% ant1_fa finds the ratio in relationship to the first and scales each
% antenna by that amplitude
ant1_fa=ant1a/ant1a*ant1;
ant2_fa=ant1a/ant2a*ant2;
ant3_fa=ant1a/ant3a*ant3;
ant4_fa=ant1a/ant4a*ant4;

% Calibrate Phases Diffrences at Broadside
% ant1_cal apply the phase diff to each of the antenna elements
% convert to radians, then multiply by e^()
ant1_cal=exp(-1i*pi/180*ant1_Phase)*ant1_fa;
ant2_cal=exp(-1i*pi/180*ant2_Phase)*ant2_fa;
ant3_cal=exp(-1i*pi/180*ant3_Phase)*ant3_fa;
ant4_cal=exp(-1i*pi/180*ant4_Phase)*ant4_fa;

% Array Vector Creation
% at every angle 1:181 we want to find the phase of each element then
% create a vector in phases in refernce to the first element
    for k=1:181
      % Find phase 
        % ant1_calmag finds the amplitude at that angle
        % ant1_FFT creates a fourier vectpr
        % ant1a is the max amplitude of that vector
        % ant1b is where it occurs
        % ant1_Phase is the phase is degrees
        ant1_calmag = max(abs(ant1_cal(k,:)));
        ant1_FFT = fft(ant1_cal(k,:));
        [ant1a,ant1b] = max(abs(ant1_FFT));
        ant1_Phase = angle(ant1_FFT(ant1b))*180/pi;

        ant2_calmag = max(abs(ant1_cal(k,:)));
        ant2_FFT = fft(ant2_cal(k,:));
        [ant2a,ant2b] = max(abs(ant2_FFT));
        ant2_Phase = angle(ant2_FFT(ant2b))*180/pi;

        ant3_calmag = max(abs(ant1_cal(k,:)));
        ant3_FFT = fft(ant3_cal(k,:));
        [ant3a,ant3b] = max(abs(ant3_FFT));
        ant3_Phase = angle(ant3_FFT(ant3b))*180/pi;

        ant4_calmag = max(abs(ant1_cal(k,:)));
        ant4_FFT = fft(ant4_cal(k,:));
        [ant4a,ant4b] = max(abs(ant4_FFT));
        ant4_Phase = angle(ant4_FFT(ant4b))*180/pi; 

        % Find the Phse Vector
        % Create the element by scaling the amplitutde and shifting it by the
        % diffrence between each element and the first element
        ARV(k,1)=ant1_calmag*exp(-1i*(ant1_Phase-ant1_Phase)*pi/180);
        ARV(k,2)=ant2_calmag*exp(-1i*(ant2_Phase-ant1_Phase)*pi/180);
        ARV(k,3)=ant3_calmag*exp(-1i*(ant3_Phase-ant1_Phase)*pi/180);
        ARV(k,4)=ant4_calmag*exp(-1i*(ant4_Phase-ant1_Phase)*pi/180);
    end
    
warning('off','all')
colorstr={[1 0 0],[255,140,0]/255,[0,0,1],[0,100,0]/255};
letterstr={'(a)','(b)','(c)','(d)'};

%Individual Instances
figure
thetaS=0; h=1; p=0;
ResponseNoInterference
legend('Location','northoutside','Orientation','horizontal','NumColumnsMode','manual','NumColumns',2)
print('NoIntResponsePolar0b','-dpdf','-bestfit')

figure
thetaS=-30; h=2; p=0;
ResponseNoInterference
thetaS=-60; h=2; p=2;
ResponseNoInterference
print('NoIntResponsePolar30_60b','-dpdf','-bestfit')

figure
thetaI1=-30;h=1;
Response1Interference
thetaI1=30;h=2;
Response1Interference
thetaI1=-60;h=3;
Response1Interference
thetaI1=60;h=4;
Response1Interference
print('OneIntResponsePolar30_60','-dpdf','-bestfit')

%%%%%%%%%%%%%%%%%%%%%%%
figure
h=1;
Response2Interference


figure
h=1;
thetaS=-10;thetaI1=-60;thetaI2=20;thetaI3=55;h=1;
Response3Interference
thetaS=-10;thetaI1=-50;thetaI2=30;thetaI3=60;h=2;
Response3Interference
thetaS=-10;thetaI1=-20;thetaI2=0;thetaI3=30;h=3;
Response3Interference
print('ThreeIntResponsePolarGoodBadAvg','-dpdf','-bestfit')

%Error Comarison
figure
spacing=1;
LMSErrorErrorOnly
GedoraErrorErrorOnly
legend('Location','north','Orientation','horizontal','NumColumnsMode','manual','NumColumns',2)
set(gca,'FontSize',14)
set(gca,'xtick',[-90 -60 -30 0 30 60 90])
print('NoIntErrorPolar','-dpdf','-bestfit')

figure
LMSError1InterferenceErrorOnly
GedoraError1InterferenceErrorOnly
legend('Location','northoutside','Orientation','horizontal','NumColumnsMode','manual','NumColumns',2)
set(gca,'FontSize',14)
set(gca,'xtick',[-90 -60 -30 0 30 60 90])
print('OneIntErrorPolarSmooth','-dpdf','-bestfit')

figure
LMSError2InterferenceErrorOnly
GedoraError2InterferenceErrorOnly
legend('Location','northoutside','Orientation','horizontal','NumColumnsMode','manual','NumColumns',2)
set(gca,'FontSize',14)
set(gca,'xtick',[-90 -60 -30 0 30 60 90])
print('TwoIntErrorPolar','-dpdf','-bestfit')