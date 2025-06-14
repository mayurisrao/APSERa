clear all;

%% use the load command to load the required mat file
load 'D:\Keerthi\PMA_tapered_circular patch'\J_top_reconv1.mat;

% load 'D:\Keerthi\Final_codes_with_sphericalsampling'\J_pc_2.mat;

 Jt=recon_J;%% look at the variable loaded and change the variable on RHS accordingly. 

nodelist=nodelist; 
freq=(2.5:.03:4)*10^9;
thetafull=-90:2:90;
phifull=0:2:180;
[t,p]=meshgrid(thetafull,phifull);
theta=thetafull;
phi=phifull;
t=t';
p=p';

ttop=(t(:));
ptop=(p(:));
r=100;  %defining the farfield radius to be 1000m
%% definition of the farfield points in cartesian 
for i=1:size(theta,2)
    for j=1:size(phi,2)
            
    xf(i,j)=r.*sind(theta(i)).*cosd(phi(j));
    yf(i,j)=r.*sind(theta(i)).*sind(phi(j));
    zf(i,j)=r.*cosd(theta(i));
    end
end
xff=xf(:);
yff=yf(:);
zff=zf(:);
w=2.*pi.*freq;
mu=4*pi*(1e-7);
epsi=8.854*1e-12;
eeta=377;
lamda0=3e8./freq;
beta0=2*pi./lamda0;
Q11t=zeros(size(xff,1),size(nodelist,1));
Q12t=zeros(size(xff,1),size(nodelist,1));
Q21t=zeros(size(xff,1),size(nodelist,1));
Q22t=zeros(size(xff,1),size(nodelist,1));

delx=nodelist(2,1)-nodelist(1,1);
dely=delx;
cellarea=delx*dely;
fieldarraytop=[xff yff zff];  

diff_matrix = reshape(fieldarraytop, [size(fieldarraytop, 1), 1, 3]) - reshape(nodelist, [1, size(nodelist, 1), 3]);
Eff=zeros(2*size(fieldarraytop,1),size(freq,2));
Efftheta=zeros(size(fieldarraytop,1),size(freq,2));
Effphi=zeros(size(fieldarraytop,1),size(freq,2));
diff_matrix2=(diff_matrix.^2);
dist = (sqrt(sum(diff_matrix2, 3))); 
dist2=dist.^2;
dist3=dist2.*dist;
dist4=dist3.*dist;
dist5=dist4.*dist;

Efftheta_=zeros(size(theta,2),size(phi,2),size(freq,2));
Effphi_=zeros(size(theta,2),size(phi,2),size(freq,2));
Emag=zeros(size(theta,2),size(phi,2),size(freq,2));
% f=45;  %% define the frequency index at which you want to calculate the farfield pattern. This should match with the frequency index at which current was calculated
for k=1:1:size(freq,2)
   
    const1=(1j*w(k)*mu./(dist3))+(3*eeta./(dist4))-(3./(1j.*w(k).*epsi.*(dist5)));
    const2=(1j*w(k)*mu./dist)+(eeta./(dist2))-(1./(1j*w(k)*epsi.*(dist3)));
    const3=exp(1j.*beta0(k).*dist);


    for i=1:size(xff,1)
        for j=1:size(nodelist,1)
            dx=diff_matrix(i,j,1);
            dy=diff_matrix(i,j,2);
            dx2=diff_matrix2(i,j,1);
            dy2=diff_matrix2(i,j,2);
            dxdy=dx*dy;
            dz=diff_matrix(i,j,3);
            convercon1=cosd(ttop(i)).*cosd(ptop(i));
            convercon2=cosd(ttop(i)).*sind(ptop(i));
            convercon3=sind(ttop(i));
            Q11t(i,j)=((convercon1*((dx2*const1(i,j))-const2(i,j)))+(convercon2*(dxdy)*const1(i,j))-(convercon3*(dx)*dz*const1(i,j)))*cellarea*const3(i,j)/4/pi;
            Q12t(i,j)=((convercon1*(dxdy)*const1(i,j))+(convercon2*((dy2*const1(i,j))-const2(i,j)))-(convercon3*(dy)*dz*const1(i,j)))*cellarea*const3(i,j)/4/pi;
            Q21t(i,j)=((cosd(ptop(i))*(dxdy)*const1(i,j))-(sind(ptop(i))*((dx2*const1(i,j))-const2(i,j))))*cellarea*const3(i,j)/4/pi;
            Q22t(i,j)=((cosd(ptop(i))*((dy2*const1(i,j))-const2(i,j)))-(sind(ptop(i))*(dxdy)*const1(i,j)))*cellarea*const3(i,j)/4/pi;


        end
    end

    H=[Q11t Q12t;Q21t Q22t];
    Eff(:,k)=H*Jt(:,k);  
    Efftheta(:,k)=Eff(1:size(fieldarraytop,1),k);
    Effphi(:,k)=Eff(size(fieldarraytop,1)+1:end,k);
    Efftheta_(:,:,k)=reshape(Efftheta(:,k),size(theta,2),size(phi,2));
    Effphi_(:,:,k)=reshape(Effphi(:,k),size(theta,2),size(phi,2));
    Emag(:,:,k) = sqrt(abs(Efftheta_(:,:,k)).^2 + abs(Effphi_(:,:,k)).^2);
    disp(k);
    
end
 %% calculation of radiation intensity

U=(Emag.^2)/377;

%% calculation of radiated power
dtheta = deg2rad(thetafull(end) - thetafull(end-1));
dphi = deg2rad(phifull(91) - phifull(90));

for  k=1:1:size(freq,2)
    Pr=0;
    for i=1:size(phifull,2)
           Pr= Pr+sum(abs(U(1:end,i,k) .* sin(deg2rad(thetafull(1:end)')) * dtheta * dphi));
    end
    Prad(k)=1*Pr;
end

%% calculation of directivity
directivity=zeros(size(thetafull,2),size(phifull,2),size(freq,2));
directivitynorm=zeros(size(thetafull,2),size(phifull,2),size(freq,2));
for i=1:1:size(freq,2)
directivity(:,:,i)=4*pi*U(:,:,i)./Prad(i);
end
 
%% normalising the directivity
for i=1:1:size(freq,2)
     for j=1:size(phifull,2)
     maxval=max(directivity(:,j,i));
     directivitynorm(:,j,i)=directivity(:,j,i)./maxval;
     end
end

    

% this is or loading the CST farfield data for comparison
allfile3=dir('D:\Keerthi\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_surfacecurrents - slots23\Export\Farfield\*[1]*.txt');
for i=1:size(allfile3,1)
    basefile3=allfile3(i).name;
    filena3=fullfile('D:\Keerthi\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_surfacecurrents - slots23\Export\Farfield\',basefile3);
    pattern = 'f=(\d+\.\d+|\d+)'; %% this is to identify the frequency values automatically from the file names
    tokens = regexp(filena3, pattern, 'tokens');
    if ~isempty(tokens)
    freq_(i) = str2double(tokens{1}{1});
    disp('Extracted frequency:');
    disp(freq_(i)); 
    end
    farf4(:,:,i)=importdata(filena3).data;
    %farf4(:,:,i)=importdata(filena3);
    %farff(:,:,i)=reshape(farf4(:,3,i),91,181);
    farff4(:,:,i)=reshape(farf4(:,3,i),181,360);
end
 %% normalising the farfield data
for k=1:1:size(freq_,2)
    for j=1:size(farff4,2)
            maxval=max(farff4(:,j,k),[],'all');
           ffnorm4(:,j,k)=farff4(:,j,k)./maxval;
    end

end


%% this is to overlay and plot multiple frequencies of CST farfield in single plot.
 figure;
for i=1:6
    % j=i;
    j=1+(i-1)*10;
    disp(freq(j));
%     polarpattern(theta,Emag(45,:,j)./max(Emag(45,:,j)));
%     polarpattern(0:180,Em4(:,135,j));
    hold on;
     polarpattern(thetafull,10.*log10(directivitynorm(:,46,j)));
end
f=44; %% change the index here to plot the below plots for a given freq


% %% this will compare and plot the MATLAB generated beam and CST generated beam
figure;
polarpattern(thetafull,10.*log10(directivitynorm(:,1,f)));
hold on;
polarpattern(0:180,10.*log10(ffnorm2(:,1,f)));
hold on;
polarpattern(thetafull,10.*log10(directivitynorm(:,46,f)));
hold on;
polarpattern(0:180,10.*log10(ffnorm2(:,91,f)));
legend('MATLAB H-plane','CST H-plane','MATLAB E-plane', 'CST E-plane');
% % 
% [tt,pp]=meshgrid(0:180,0:359);
% tt=tt';
% pp=pp';
% %%rms error between the 3d beams
% for  k=1:1:size(freq,2)
% for i=1:size(phi,2)
%     di(:,i,k)=interp1(theta(46:end),directivitynorm(46:end,i,k),0:90,'cubic');
%     di2(:,i,k)=interp1(theta(46:end),flip(directivitynorm(1:46,i,k)),0:90,'cubic');
% end
% 
% di_f(:,:,k)=[di(:,:,k) di2(:,2:end,k)];
% end
% for  k=1:1:size(freq,2)
% for i=1:size(di_f,1)
%     di_final(i,:,k)=interp1(0:2:360,di_f(i,:,k),0:359,'cubic');
% end
% rmserror(k)=rms(di_final(:,:,k)-ffnorm(1:91,:,k),[1 2]);
% end 

[Gx,Gy,Gz]=gradient(ffnorm4(:,:,:),1,1,1);
chrom=chromaticity1(Gz(:,:,1:end)); 
 
for i=1:size(freq_,2)-1
    [ssimval,ssimmap]=ssim(ffnorm4(:,:,i).*1e3,ffnorm4(:,:,i+1).*1e3,'Radius',1);   %% the directivity values are scaled by 1e3 because the similarity index recognizes the differences more prominently
     simval7(i)=ssimval;
     ssmap(:,:,i)=ssimmap;
end
% disp('The redisual error between the CST beam and MATLAB beam is');
% disp(rmserror);
% Emaghi=Emag(:,:,5:14);
% directivitynormhi=directivitynorm(:,:,5:14);
% Effphihi_=Effphi_(:,:,5:14);
% Effthetahi_=Efftheta_(:,:,5:14);
% rmserrorhi=rmserror(5:14);
%% save the farfield data (Emag, Etheta, Ephi, directivity) in a .mat file
        save('D:\Keerthi\PMA_tapered_circular patch\patch_recontopffv1.mat','Emag','directivitynorm','Effphi_','Efftheta_','nodelist');
%      save('D:\Keerthi\Final_codes_with_sphericalsampling\farield_10_20.mat','ffnorm');
                      % load('D:\Keerthi\Final_codes_with_sphericalsampling\farfield_pc2.mat');
%    % r  
%    % ms2=rmserror;
%  rms2=rmserror;
% directivitynorm(:,:,10:31)=directivitynormhi;
% Emag(:,:,10:31)=Emaghi;
% Effphi_(:,:,10:31)=Effphihi_;
% Efftheta_(:,:,10:31)=Effthetahi_;
% rmserror(10:31)=rms2(10:31);


 load('D:\Keerthi\PMA_tapered_circular patch\patch_recontopff.mat');
Emagtop=Emag;
 load ('D:\Keerthi\PMA_tapered_circular patch\patch_bottomff.mat');
for i=1:size(freq,2)
Emagbottom(:,:,i)=Emag(:,:,i)./max(Emag(:,:,i),[],'all');
end
thetafull=0:2:180;
phifull=0:2:360;
% Emag=[[Emagbottom(1:46,:,:);Emagtop(2:46,:,:)] Emagbottom(47:end,:,:)];
Emag1=[Emagtop(46:end,:,:); Emagbottom(48:end,:,:)];
Emag2=[flip(Emagtop(2:46,2:end,:),1); flip(Emagbottom(1:46,2:end,:),1)];
Emag=[Emag1 Emag2];
% 
[tt,pp]=meshgrid(thetafull,phifull);
tt=tt';
pp=pp';
tfull=tt(:);
pfull=pp(:);
directory_path = 'D:\Keerthi\PMA_tapered_circular patch\farff'; 
ffarray=zeros(size(tfull,1),3);
for i=1:size(freq,2)

    temp=directivitynorm(:,:,i);

    dfull(:,i)=temp(:);
    ffarray=[tfull pfull dfull(:,i)];
    filename = sprintf('farfield (f=%.3f).txt', freq(i).*1e-9); 
    full_file_path = fullfile(directory_path, filename);
    fid = fopen(full_file_path, 'w');  % 'w' for writing (creates or overwrites the file)
    % Write the matrix to the file, separating columns by tabs
    for i = 1:size(ffarray, 1)  % Loop through each row
    fprintf(fid, '%d\t%d\t%f\n', ffarray(i,:));  % Print each row with tabs separating columns
    end

% Close the file
    fclose(fid);
end
% 
% % Emag=[Emagbottom(1:45,:,:);Emagtop(1:end,:,:);Emagbottom(48:end,:,:)];


a=importdata('C:\Users\APSERa\Documents\vna_meas_data\5may2025_cryosw15.csv');
s=a.data;
s11=s(:,2)+1j.*s(:,3);
s21=s(:,4)+1j.*s(:,5);
smag1=20.*log10(abs(s11));
smag2=20.*log10(abs(s21));
freq=cell2array(a.textdata);
figure;
plot(smag1);
hold on;
plot(smag2);
