clear all;
freq=2.605:0.015:4;
%% simulated data import
for i=1:size(freq,2)
    basefile1=sprintf('farfield (f=%0.3f) [1].txt',freq(i));
    filena1=fullfile('F:\rri files\APSERA\anechoic chamber meas\with sma',basefile1);
    c2(:,:,i)=importdata(filena1).data;
    for ii=0:359
    jj=(ii*180)+(ii+1);
    simulgain(:,ii+1,i)=c2(jj:jj+180,3,i);
    end
end
%% for replicating the theta value from 0 to 180 to -180 to 0
 simulgain1=flip(simulgain,1);
 simulgaintotal=[simulgain1 ;simulgain(2:end,:,:)]; 


 for i=1:size(freq,2)
    basefile1=sprintf('farfield (f=%0.3f) [2].txt',freq(i));
    filena1=fullfile('F:\rri files\APSERA\anechoic chamber meas\with sma',basefile1);
    c3(:,:,i)=importdata(filena1).data;
    for ii=0:359
    jj=(ii*180)+(ii+1);
    simulgain_ant2(:,ii+1,i)=c3(jj:jj+180,3,i);
    end
end
%% for replicating the theta value from 0 to 180 to -180 to 0
 simulgain_ant21=flip(simulgain_ant2,1);
 simulgaintotal_ant2=[simulgain_ant21 ;simulgain_ant2(2:end,:,:)]; 

%  polarpattern(-180:180,simulgaintotal(:,1,80));
%  hold on
%   polarpattern(-180:180,simulgaintotal_ant2(:,1,80));

%%% measured data antenna 1 copol
 eplane1=importdata("F:\rri files\APSERA\anechoic chamber meas\config1\antenna 1\eplane copol.txt");
 hplane1=importdata("F:\rri files\APSERA\anechoic chamber meas\config1\antenna 1\hplane copol.txt");
 eplanecross1=importdata("F:\rri files\APSERA\anechoic chamber meas\config1\antenna 1\eplane crosspol.txt");
 hplanecross1=importdata("F:\rri files\APSERA\anechoic chamber meas\config1\antenna 1\hplane crosspol.txt");
 eplane2=importdata("F:\rri files\APSERA\anechoic chamber meas\config1\antenna 2\eplane copol.txt");
 hplane2=importdata("F:\rri files\APSERA\anechoic chamber meas\config1\antenna 2\hplane copol.txt");
 eplanecross2=importdata("F:\rri files\APSERA\anechoic chamber meas\config1\antenna 2\eplane crosspol.txt");
 hplanecross2=importdata("F:\rri files\APSERA\anechoic chamber meas\config1\antenna 2\hplane crosspol.txt");
 horn=importdata("F:\rri files\APSERA\anechoic chamber meas\horndata.csv");
 hornhorn=importdata("F:\rri files\APSERA\anechoic chamber meas\Horn_horn.txt");
% 
% 
% 
% 
% hornhorn=aa.data;
hornf=horn(:,1);
gainhorn=horn(:,2);
freq_=2:0.015:5;
gain_interp=interp1(hornf,gainhorn,freq_);
plot(freq_,gain_interp);

%% antenna 1 meas
a=eplane1.data;
theta=1:361;
freq1=a(:,2);
gain=a(:,3);

b=hplane1.data;
theta1=b(:,1);
gain1=b(:,3);
for j=1:361
    
%     for i=1:200
        gain_ep(j,:)=gain(j+(200*(j-1)):j+(200*(j)),1);
        gain_hp(j,:)=gain1(j+(200*(j-1)):j+(200*(j)),1);
        
%     end
end
gainantep=gain_ep-repmat(hornhorn(:,2)',361,1)+repmat(gain_interp,361,1);
gainanthp=gain_hp-repmat(hornhorn(:,2)',361,1)+repmat(gain_interp,361,1);

% antenna 1 cross meas
ac=eplanecross1.data;
% theta=a(:,1);
% freq1=a(:,2);
gaincross1=ac(:,3);

bc=hplanecross1.data;
% theta1=b(:,1);
gaincross1_=bc(:,3);
for j=1:361
%     for i=1:200
        gain_epcross(j,:)=gaincross1(j+(200*(j-1)):j+(200*(j)),1);
        gain_hpcross(j,:)=gaincross1_(j+(200*(j-1)):j+(200*(j)),1);
        
%     end
end
gainantepcross1=gain_epcross-repmat(hornhorn(:,2)',361,1)+repmat(gain_interp,361,1);
gainanthpcross1=gain_hpcross-repmat(hornhorn(:,2)',361,1)+repmat(gain_interp,361,1);

%% antenna 2 meas
a_=eplane2.data;
theta=a_(:,1);
% freq1=a_(:,2);
gain2=a_(:,3);

b_=hplane2.data;
% theta1=b(:,1);
gain3=b_(:,3);
for j=1:361
%     for i=1:200
        gain_ep2(j,:)=gain2(j+(200*(j-1)):j+(200*(j)),1);
        gain_hp2(j,:)=gain3(j+(200*(j-1)):j+(200*(j)),1);
        
%     end
end
gainantep2=gain_ep2-repmat(hornhorn(:,2)',361,1)+repmat(gain_interp,361,1);
gainanthp2=gain_hp2-repmat(hornhorn(:,2)',361,1)+repmat(gain_interp,361,1);

%% antenna 2 cross meas
ac2=eplanecross2.data;
% theta=a_(:,1);
% freq1=a_(:,2);
gaincross2=ac2(:,3);

bc2=hplanecross2.data;
% theta1=b(:,1);
gaincross2_=bc2(:,3);
for j=1:361
%     for i=1:200
        gain_epcross2(j,:)=gaincross2(j+(200*(j-1)):j+(200*(j)),1);
        gain_hpcross2(j,:)=gaincross2_(j+(200*(j-1)):j+(200*(j)),1);
        
%     end
end
gainantepcross2=gain_epcross2-repmat(hornhorn(:,2)',361,1)+repmat(gain_interp,361,1);
gainanthpcross2=gain_hpcross2-repmat(hornhorn(:,2)',361,1)+repmat(gain_interp,361,1);



%% eplane in simulation is phi=135 degree for antenna 1
%% hplane in simulation is phi=225 or 45 degree for antenna 1
%% 2.608GHz is in 8th index in simulation
%% 2.8GHz is in 21st index in simulation
%%3.01GHz is in 31st index
%%3.5GHz is in 68th index
%% 4GHz is in 101th index
plotf=[8 21 54 68 82];
%% extracting the freq data from simulation
for ii=1:size(plotf,2)
normant1ep(:,ii)=simulgaintotal(:,135,plotf(ii));
normant1hp(:,ii)=simulgaintotal(:,225,plotf(ii));
normant2ep(:,ii)=simulgaintotal_ant2(:,45,plotf(ii));
normant2hp(:,ii)=simulgaintotal_ant2(:,135,plotf(ii));

end

%% max normalising the data
normant1ep_max=(normant1ep-repmat(max(normant1ep,[],1),361,1));
normant1hp_max=(normant1hp-repmat(max(normant1hp,[],1),361,1));
normant2ep_max=(normant2ep-repmat(max(normant2ep,[],1),361,1));
normant2hp_max=(normant2hp-repmat(max(normant2hp,[],1),361,1));
gainant1ep_max=(gainantep-repmat(max(gainantep,[],1),361,1));
gainant1hp_max=(gainanthp-repmat(max(gainanthp,[],1),361,1));
gainant1epcross_max=(gainantepcross1-repmat(max(gainantep,[],1),361,1));
gainant1hpcross_max=(gainanthpcross1-repmat(max(gainanthp,[],1),361,1));
gainant2ep_max=(gainantep2-repmat(max(gainantep2,[],1),361,1));
gainant2hp_max=(gainanthp2-repmat(max(gainanthp2,[],1),361,1));
gainant2epcross_max=(gainantepcross2-repmat(max(gainantep2,[],1),361,1));
gainant2hpcross_max=(gainanthpcross2-repmat(max(gainanthp2,[],1),361,1));

%% plotting for co and cross of measured
figure;
polarpattern(-180:180,normant1ep_max(:,1),-180:180,gainant1ep_max(:,41));
figure;
polarpattern(-180:180,gainant1hp_max(:,55),-180:180,gainant2hp_max(:,55));
figure;
polarpattern(-180:180,gainant1hp_max(:,68),-180:180,gainant2hp_max(:,68));
figure;
polarpattern(-180:180,gainant1hp_max(:,101),-180:180,gainant2hp_max(:,101));
figure;
polarpattern(-180:180,gainant1hp_max(:,128),-180:180,gainant2hp_max(:,128));
figure;
polarpattern(-180:180,gainant2hp_max(:,50),-180:180,gainant2hp_max(:,55),-180:180,gainant2hp_max(:,68),-180:180,gainant2hp_max(:,101),-180:180,gainant2hp_max(:,128));

%%PATTERNFROMSLICES function takes theta as input, but the pattern should
%%have maximum towards x axis, so the 'thet' is fed as elevation angles...and
%%phi varies from -180 to 180. The function also takes theta values greater
%%than 180 only. The crossweighted method is more accurate than the summing
%%method.
thet=90-(-180:179);
for kkk=1:101
abc(:,:,kkk)=patternFromSlices(gainant2hp_max(1:360,kkk+35),thet,gainant2ep_max(1:360,kkk+35),-180:179,'method','CrossWeighted');
rotpattern(:,:,kkk)=rotpat(abc(:,:,kkk)',-180:179,-90:90,rotx(90)*rotz(90),-inf); %% Here the pattern maxima from x axis has to be rotated to z axis
rotapattern1(:,:,kkk)=flip(rotpattern(:,:,kkk),1); %% the rotatedpattern is stored such a way that the maxima occurs at theta=180, so it is flipped here to make maxima at theta=0
a1=rotapattern1(:,:,kkk); %% some values are not reconsuctred properly and they have been assigned -inf. These values are then replaced by -40dB.
a1(a1==-inf)=-40;

%% the next portion is to clean out some improperly reconstructed values in theta=90 plane. These values are identified to be around nulls. 
rotapattern1(:,:,kkk)=a1;
rotapattern1(:,:,kkk)=(rotapattern1(:,:,kkk)-repmat(max(rotapattern1(:,:,kkk),[],1),181,1));
ai=gradient(rotapattern1(90,:,kkk));
    [m1,ind]=max(abs(ai(1:180)));
    [m,ind1]=max(abs(ai(180:end)));
rotapattern1(90,ind-3:ind+3,kkk)=min(rotapattern1(90,:,kkk)); %% the outliers are assigned the minimum value of array
rotapattern1(90,ind1-3+180:ind1+3+180,kkk)=min(rotapattern1(90,:,kkk)); %% the outliers are assigned minimum values of array
end

%% writing the values into a proper format compatible with python pipeline
for kkk=1:101
kk=1;
for ii=1:360
        abc_pip(kk:kk+180,3,kkk)=rotapattern1(:,ii,kkk);
        abc_pip(kk:kk+180,1,kkk)=0:180;
        abc_pip(kk:kk+180,2,kkk)=ii-1;
        kk=kk+1+180;
end
end

%% writing the reconstructed 3d pattern to frequency files ( Apsera pipeline compatible format)
for i=1:101
j=i+35;
patternfile=sprintf('%0.3fghz.txt',freq_(j));
patterna=fullfile('F:\rri files\APSERA\anechoic chamber meas\config1\3d pattern\antenna2',patternfile);
fileID = fopen(patterna,'w');
fprintf(fileID,'%d %d %0.4f \r\n',abc_pip(:,:,i)');
end


%% Sparameter plots
cd=importdata('F:\rri files\APSERA\anechoic chamber meas\keerthi s11_.txt');
dd=importdata('F:\rri files\APSERA\anechoic chamber meas\keeerthi s22_.txt');
s11plot=smooth(cd(:,2),50);
figure(8);
plot(freq_,smooth(cd(:,2),50),freq_,smooth(dd(:,2),50));
axis([2.5 4 -20 0]);

ee=importdata('F:\rri files\APSERA\anechoic chamber meas\s2p.csv');
fr=ee(:,1);
s11_sim=ee(:,2);
s21_sim=ee(:,3);
s22_sim=ee(:,4);
plot(fr,s11_sim,freq_,smooth(cd(:,2),50),fr,s22_sim,freq_,smooth(dd(:,2),50));
axis([2.5 4 -20 0]);
