%##########################################################################
%                     Wavelet-epileptic-seizure-EEG                       #
%                      Presented by: Reza Saadatyar                       #
%                        File Name: chb01_26.edf                          #
%                        File Start Time: 12:34:22                        #
%                        File End Time: 13:13:07                          #
%                        Seizure Start Time: 1862 seconds                 #
%                        Seizure End Time: 1963 seconds                   #
%                              Fs = 256 Hz                                #
%##########################################################################
% Please Run Code Main.m, the Software will run automatically and there is no need to run others code
clc;clear;close all;
%% Load data
Fs = 256;
data=load('data.mat');
data=data.data;
data=data(1,440250:490761);Time= (1:length(data))/Fs;
dataNonseizures=data(:,1:25250); Time1= (1:length(dataNonseizures))/Fs;% non -seizures
dataseizures=data(:,25250:50499); % seizures
figure(1);
subplot(2,2,[1 2]);
plot(Time,data); 
xlabel('Time(Sec)','FontSize',12,'FontWeight','bold','FontName','Times New Roman'); 
ylabel('Voltage(\muV)','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
title('Epilepsy EEG(seizures & Non-seizures); Sampling Rate: 256 Hz','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
legend({'Epilepsy EEG(seizures & Non-seizures)'},'FontSize',9,'FontWeight','bold','FontName','Times New Roman');
grid on;grid minor;
subplot(223)
plot(Time1,dataNonseizures); 
xlabel('Time(Sec)','FontSize',12,'FontWeight','bold','FontName','Times New Roman'); 
ylabel('Voltage(\muV)','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
legend({'Epilepsy EEG(Non-seizures)'},'FontSize',9,'FontWeight','bold','FontName','Times New Roman');
grid on;grid minor;
subplot(224)
plot(Time1,dataseizures); 
xlabel('Time(Sec)','FontSize',12,'FontWeight','bold','FontName','Times New Roman'); 
ylabel('Voltage(\muV)','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
legend({'Epilepsy EEG(seizures)'},'FontSize',9,'FontWeight','bold','FontName','Times New Roman');
grid on;grid minor;
%% design wavelet
wname='db4';
nLevel=4;
[a,d]=swt(data,nLevel,wname);
figure(2);
subplot(nLevel+1,2,[1 2]);
plot(Time,data(1,:)); ylabel('Epilepsy EEG');
ylabel('\muV','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
title('Epilepsy EEG(seizures & Non-seizures)','FontName','Times New Roman','FontSize',12,'FontWeight','bold');
grid on;grid minor;
c=2;
for i=1:nLevel
    c=c+1;
    subplot(nLevel+1,2,c);
    plot(Time,a(i,:),'c');
    ylabel(['a_{' num2str(i) '}'],'FontSize',12,'FontWeight','bold','FontName','Times New Roman');
    grid on;grid minor;
     if c==9
       xlabel('Time(Sec)','FontSize',12,'FontWeight','bold','FontName','Times New Roman'); 
    end
    c=c+1;
    subplot(nLevel+1,2,c);
    plot(Time,d(i,:),'r');
    ylabel(['d_{ ' num2str(i) '}'],'FontSize',12,'FontWeight','bold','FontName','Times New Roman');
    grid on;grid minor;
    if c==10
    xlabel('Time(Sec)','FontSize',12,'FontWeight','bold','FontName','Times New Roman'); 
    end
end
%% coffiction wavelet
[c,l] = wavedec(data,nLevel,wname);
ca4=appcoef(c,l,wname,4);
[cd1,cd2,cd3,cd4]=detcoef(c,l,[1,2,3,4]);
%% Threshould
alpha=3;
Y1=find(abs(cd1)<= alpha);cd1(Y1)=0;
Y2=find(abs(cd2)<= alpha);cd2(Y2)=0;
Y3=find(abs(cd3)<= alpha);cd3(Y3)=0;
Y4=find(abs(cd4)<= alpha);cd4(Y4)=0;
X4=find(abs(ca4)<=alpha);ca4(X4)=0;
yy=cat(2,ca4,cd4,cd3,cd2,cd1);
R1= waverec(yy,l,wname);
%% Round Coefficient
Y11=round(cd1);
Y22=round(cd2);
Y33=round(cd3);
Y44=round(cd4);
X44=round(ca4);
yy=cat(2,X44,Y44,Y33,Y22,Y11);
R11= waverec(yy,l,wname);
%% Quantization 
partition =[min(R11):70:max(R11)];%#ok
codebook =[min(R11)-50:70:max(R11)];%#ok
[partition2,codebook2]=lloyds(R11,codebook);%optimize, using codebook 
[index,quants,distor]=quantiz(R11,partition,codebook);
[index2,quant2,distor2]=quantiz(R11,partition2,codebook2);
% figure();
% subplot(211);
% plot(quant2,'r-');
% xlabel('Time(Sec)','FontSize',12,'FontWeight','bold','FontName','Times New Roman'); 
% ylabel('Voltage(\muV)','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
% title('Quantizing Recognition Epilepsy EEG Signal','FontName','Times New Roman','FontSize',16,'FontWeight','bold');
% set(gca,'xticklabel',{0:40:240});
% grid on;grid minor;
% subplot(212);
% plot(quant2(1,1:500),'g-');
% xlabel('Time(msec)','FontSize',12,'FontWeight','bold','FontName','Times New Roman'); 
% ylabel('Voltage(\muV)','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
% title('Quantizing Recognition Epilepsy EEG Signal','FontName','Times New Roman','FontSize',16,'FontWeight','bold');
% set(gca,'xticklabel',{0:100:1000});
% grid on;grid minor;
%% Arithmetic coding
quant2=quant2/800;
quant2=quant2+1+abs(min(quant2));
quant2=round(quant2);
counts = [50 26 12 12];
code = arithenco(quant2,counts);
SizeOrginalSignal=size(R11,2);
SizeComperseSignal = size(code,2);
CR=SizeComperseSignal/SizeOrginalSignal;
R=(1/CR)*100;
R=100-R;

figure(3);
plot(Time,data,'r');
ylabel('\muV','FontSize',12,'FontWeight','bold','FontName','Times New Roman');
grid on;grid minor;hold on
plot(Time,R1,'m');
Energy=(norm(R1,2)/norm(data,2))*100;% Energy
xlabel('Time(Sec)','FontSize',12,'FontWeight','bold','FontName','Times New Roman'); 
title([' Recognition Epilepsy EEG Signal with Energy :',num2str(Energy),'; CR:',num2str(CR),', R:%',num2str(R),';  \alpha=2.5'],'FontName','Times New Roman','FontSize',16,'FontWeight','bold');
legend('Raw EEG','Recognition EEG')
%% feature standardized
dataseizures=R11(:,1:25250);
dataNonseizures=R11(:,25250:50499);
k=0;
for i=1:12625:25250
    k=k+1;
    EEGdataseizures(:,k)=dataseizures(i:i+12624);%#ok
    EEGdataNonseizures(:,k)=dataNonseizures(i:i+12624);%#ok
end
X=[EEGdataseizures-350;EEGdataseizures+100];
N=size(EEGdataseizures,1);
y=[ones(N,1);2*ones(N,1)];
MeanEEGclass1=mean(EEGdataseizures);% non -seizures
StdEEGclass1=std(EEGdataseizures);
EEGclass1=((EEGdataseizures-MeanEEGclass1)./StdEEGclass1);
MeanEEGclass2=mean(EEGdataNonseizures);% non -seizures
StdEEGclass2=std(EEGdataNonseizures);
EEGclass2=((EEGdataNonseizures-MeanEEGclass2)./StdEEGclass2);
%% Classification
KFold=5;
%% SVM
SVM = templateSVM('Standardize',1,'KernelFunction','gaussian');
indices = crossvalind('Kfold',y,KFold);Perfomance=[];
for i = 1:KFold
    test =(indices == i);train1 = ~test;
    TrainInputs=X(train1,:);TrainTargets=y(train1,:);
    TestInputs=X(test,:);labelTargetTest=y(test,:);
    Mdl = fitcecoc(TrainInputs,TrainTargets,'Learners',SVM,'FitPosterior',1,...
        'ClassNames',1:2,'Verbose',2);label = predict(Mdl,TestInputs);
    [Acc,Sen,Spe]=ConMax(labelTargetTest,label);Perfomance=[Perfomance;Acc Sen Spe]; %#ok   
end
clc;Perfomance=mean(Perfomance,1);disp('Svm: ACC SEN SPE')
disp(Perfomance)
%% KNN
Perfomance=[];indices = crossvalind('Kfold',y,KFold); 
for i = 1:KFold
   test =(indices == i);train1 = ~test;
    TrainInputs=X(train1,:);TrainTargets=y(train1,:);
    TestInputs=X(test,:);labelTargetTest=y(test,:);
    Mdl = fitcknn(TrainInputs,TrainTargets,'NumNeighbors',5,'Standardize',1);
    label = predict(Mdl,TestInputs);
    [Acc,Sen,Spe]=ConMax(labelTargetTest,label);
    Perfomance=[Perfomance;Acc Sen Spe]; %#ok
end
Perfomance=mean(Perfomance,1);disp('KNN: ACC SEN SPE');disp(Perfomance)
%% Naive Bayesian
Perfomance=[];indices = crossvalind('Kfold',y,KFold);
for i = 1:KFold
    test =(indices == i);train1 = ~test;
    TrainInputs=X(train1,:);TrainTargets=y(train1,:);
    TestInputs=X(test,:);labelTargetTest=y(test,:);
    CMdl=fitcnb(TrainInputs,TrainTargets);
    label = predict(CMdl,TestInputs);[Acc,Sen,Spe]=ConMax(labelTargetTest,label);
    Perfomance=[Perfomance;Acc Sen Spe]; %#ok
end
Perfomance=mean(Perfomance,1);disp('Naive Bayesian: ACC SEN SPE');disp(Perfomance)
%% MLP
Label=ind2vec(y');
Perfomance=[];indices = crossvalind('Kfold',vec2ind(Label),KFold); 
for i = 1:KFold
    test =(indices == i);train1 = ~test;
    TrainInputs=X(train1,:)';TrainTargets=Label(:,train1);
    TestInputs=X(test,:)';labelTargetTest=Label(:,test);
    net = feedforwardnet(3,'trainlm');
    c=net.numLayers;net.layers{c}.transferFcn ='purelin';
    net.trainParam.epochs=100;
     net=train(net,TrainInputs,full(TrainTargets));
    label = net(TestInputs);
    [Acc,Sen,Spe]=ConMax(vec2ind(labelTargetTest),vec2ind(label));
    Perfomance=[Perfomance;Acc Sen Spe]; %#ok
end
Perfomance=mean(Perfomance,1);disp('MLP: ACC SEN SPE');disp(Perfomance)





