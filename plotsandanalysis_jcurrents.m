
%%% This initial portion of the code takes the top and bottom currents,
%%% their modified versions, reconstructed versions and find the
%%% corresponding fft, after reshaping.  The code is flexible 


%% This portion is for the top plane currents
 load('F:\rri files\APSERA\surface current data\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_140mm\Final_codes_with_sphericalsampling\patch_original\current_patch.mat');
 Jt=J;
 load('F:\rri files\APSERA\surface current data\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_140mm\Final_codes_with_sphericalsampling\Patch_modified\current_patch.mat');
  Jchange=J;
 nodelist=nodelist;
load('F:\rri files\APSERA\surface current data\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_140mm\Final_codes_with_sphericalsampling\patch_original\J_top_reconv1.mat');
freq=2.5:0.03:4;
delf=freq(2)-freq(1);
% This loop is for the top plane original and reconstructed currents
for i=1:size(freq,2)
    Jfinal=Jt;
    Jx1(:,:,i)=reshape(Jfinal(1:size(nodelist,1),i),sqrt(size(nodelist,1)),sqrt(size(nodelist,1)));
    Jx1(:,:,i)=transpose(Jx1(:,:,i));
    Jy1(:,:,i)=reshape(Jfinal(size(nodelist,1)+1:end,i),sqrt(size(nodelist,1)),sqrt(size(nodelist,1)));
    Jy1(:,:,i)=transpose(Jy1(:,:,i));
    Jmag(:,:,i)=sqrt(abs(Jx1(:,:,i)).^2+abs(Jy1(:,:,i)).^2);
    Jangle(:,:,i)=atan2(real(Jy1(:,:,i)),real(Jx1(:,:,i)));
    % Jangle(:,:,i)=angle(Jx1(:,:,i)+Jy1(:,:,i));
    Jcomplex(:,:,i)=Jmag(:,:,i).*exp(1j.*Jangle(:,:,i));

    Jfinal1=recon_J;
    Jx2(:,:,i)=reshape(Jfinal1(1:size(nodelist,1),i),sqrt(size(nodelist,1)),sqrt(size(nodelist,1)));
    Jx2(:,:,i)=transpose(Jx2(:,:,i));
    Jy2(:,:,i)=reshape(Jfinal1(size(nodelist,1)+1:end,i),sqrt(size(nodelist,1)),sqrt(size(nodelist,1)));
    Jy2(:,:,i)=transpose(Jy2(:,:,i));
    % Jmag2(:,:,i)=abs(Jx2(:,:,i)+Jy2(:,:,i));
    Jmag2(:,:,i)=sqrt(abs(Jx2(:,:,i)).^2+abs(Jy2(:,:,i)).^2);
    Jangle2(:,:,i)=atan2(real(Jy2(:,:,i)),real(Jx2(:,:,i)));
    % Jangle2(:,:,i)=angle(Jx2(:,:,i)+Jy2(:,:,i));
    Jcomplex2(:,:,i)=Jmag2(:,:,i).*exp(1j.*Jangle2(:,:,i));
    
end

% This portion is for the top plane modified design currents.
for i=1:size(Jchange,2)
    Jfinal2=Jchange;
    Jx3(:,:,i)=reshape(Jfinal2(1:size(nodelist,1),i),sqrt(size(nodelist,1)),sqrt(size(nodelist,1)));
    Jx3(:,:,i)=transpose(Jx3(:,:,i));
    Jy3(:,:,i)=reshape(Jfinal2(size(nodelist,1)+1:end,i),sqrt(size(nodelist,1)),sqrt(size(nodelist,1)));
    Jy3(:,:,i)=transpose(Jy3(:,:,i));
    Jmag3(:,:,i)=sqrt(abs(Jx3(:,:,i)).^2+abs(Jy3(:,:,i)).^2);
    Jangle3(:,:,i)=atan2(real(Jy3(:,:,i)),real(Jx3(:,:,i)));
    % Jangle(:,:,i)=angle(Jx1(:,:,i)+Jy1(:,:,i));
    Jcomplex3(:,:,i)=Jmag3(:,:,i).*exp(1j.*Jangle3(:,:,i));

       
end

%% for plotting the currents (Here we are reshaping the nodes to plot them as images. Note that we are taking only the range of values between the -0.7 to 0.7 ranges
nodelistt=reshape(nodelist,sqrt(size(nodelist,1)),sqrt(size(nodelist,1)),3);
nodelistt(:,:,1)=transpose(nodelistt(:,:,1));
nodelistt(:,:,2)=transpose(nodelistt(:,:,2));
nodelistt_=nodelistt(15:80-15,15:80-15,:);
nodelisthi3t=reshape(nodelistt_,size(nodelistt_,1)*size(nodelistt_,2),3);

%% This portion is for performing FFT on all the top plane data
fftsizetop=512;

for i=1:size(freq,2)

Jxfftop(:,:,i)=(fftshift(fft2(Jx1(:,:,i),fftsizetop,fftsizetop))); %% this is 2D spatial FFT
Jyfftop(:,:,i)=(fftshift(fft2(Jy1(:,:,i),fftsizetop,fftsizetop))); 


Jxfftop2(:,:,i)=(fftshift(fft2(Jx2(:,:,i),fftsizetop,fftsizetop))); %% this is 2D spatial FFT
Jyfftop2(:,:,i)=(fftshift(fft2(Jy2(:,:,i),fftsizetop,fftsizetop))); 


Jmagfft(:,:,i)=sqrt(abs(Jxfftop(:,:,i)).^2+abs(Jyfftop(:,:,i)).^2);
Jmagfft2(:,:,i)=sqrt(abs(Jxfftop2(:,:,i)).^2+abs(Jyfftop2(:,:,i)).^2);


power2d1(:,:,i)=(Jmagfft(:,:,i).^2)./max(Jmagfft(:,:,i).^2,[],'all');
power2d2(:,:,i)=(Jmagfft2(:,:,i).^2)./max(Jmagfft2(:,:,i).^2,[],'all');


end

for i=1:size(Jchange,2)
    Jxfftop3(:,:,i)=(fftshift(fft2(Jx3(:,:,i),fftsizetop,fftsizetop))); %% this is 2D spatial FFT
    Jyfftop3(:,:,i)=(fftshift(fft2(Jy3(:,:,i),fftsizetop,fftsizetop))); 
    Jmagfft3(:,:,i)=sqrt(abs(Jxfftop3(:,:,i)).^2+abs(Jyfftop3(:,:,i)).^2);
    power2d3(:,:,i)=(Jmagfft3(:,:,i).^2)./max(Jmagfft3(:,:,i).^2,[],'all');
end


%% This portion is for the bottom plane currents original, reconstructed, and modified
load('F:\rri files\APSERA\surface current data\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_140mm\Final_codes_with_sphericalsampling\original_design_full3d\current.mat');
Jbottom=Jhi3;
 load('F:\rri files\APSERA\surface current data\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_140mm\slots22\current_bottom.mat');
  Jchangebottom=Jhi3;
 nodelist=nodelisthi3;
load('F:\rri files\APSERA\surface current data\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_140mm\Final_codes_with_sphericalsampling\original_design_full3d\J_bottompc.mat');
recon_Jbottom=recon_Jbottom;

delf=freq(2)-freq(1);
for i=1:size(freq,2)
    Jfinalb=Jbottom;
    Jxb1(:,:,i)=reshape(Jfinalb(1:size(nodelist,1),i),sqrt(size(nodelist,1)),sqrt(size(nodelist,1)));
    Jxb1(:,:,i)=transpose(Jxb1(:,:,i));
    Jyb1(:,:,i)=reshape(Jfinalb(size(nodelist,1)+1:end,i),sqrt(size(nodelist,1)),sqrt(size(nodelist,1)));
    Jyb1(:,:,i)=transpose(Jyb1(:,:,i));
    Jmagb(:,:,i)=sqrt(abs(Jxb1(:,:,i)).^2+abs(Jyb1(:,:,i)).^2);
    Jangleb(:,:,i)=atan2(real(Jyb1(:,:,i)),real(Jxb1(:,:,i)));
    % Jangle(:,:,i)=angle(Jx1(:,:,i)+Jy1(:,:,i));
    Jcomplexb(:,:,i)=Jmagb(:,:,i).*exp(1j.*Jangleb(:,:,i));

    Jfinalb1=recon_J;
    Jxb2(:,:,i)=reshape(Jfinalb1(1:size(nodelist,1),i),sqrt(size(nodelist,1)),sqrt(size(nodelist,1)));
    Jxb2(:,:,i)=transpose(Jxb2(:,:,i));
    Jyb2(:,:,i)=reshape(Jfinalb1(size(nodelist,1)+1:end,i),sqrt(size(nodelist,1)),sqrt(size(nodelist,1)));
    Jyb2(:,:,i)=transpose(Jyb2(:,:,i));
    % Jmag2(:,:,i)=abs(Jx2(:,:,i)+Jy2(:,:,i));
    Jmagb2(:,:,i)=sqrt(abs(Jxb2(:,:,i)).^2+abs(Jyb2(:,:,i)).^2);
    Jangleb2(:,:,i)=atan2(real(Jyb2(:,:,i)),real(Jxb2(:,:,i)));
    % Jangle2(:,:,i)=angle(Jx2(:,:,i)+Jy2(:,:,i));
    Jcomplexb2(:,:,i)=Jmagb2(:,:,i).*exp(1j.*Jangleb2(:,:,i));
    
end
for i=1:size(Jchangebottom,2)
    Jfinalb2=Jchangebottom;
    Jxb3(:,:,i)=reshape(Jfinalb2(1:size(nodelist,1),i),sqrt(size(nodelist,1)),sqrt(size(nodelist,1)));
    Jxb3(:,:,i)=transpose(Jxb3(:,:,i));
    Jyb3(:,:,i)=reshape(Jfinalb2(size(nodelist,1)+1:end,i),sqrt(size(nodelist,1)),sqrt(size(nodelist,1)));
    Jyb3(:,:,i)=transpose(Jyb3(:,:,i));
    Jmagb3(:,:,i)=sqrt(abs(Jxb3(:,:,i)).^2+abs(Jyb3(:,:,i)).^2);
    Jangleb3(:,:,i)=atan2(real(Jyb3(:,:,i)),real(Jxb3(:,:,i)));
    % Jangle(:,:,i)=angle(Jx1(:,:,i)+Jy1(:,:,i));
    Jcomplexb3(:,:,i)=Jmagb3(:,:,i).*exp(1j.*Jangleb3(:,:,i));

       
end



%% This is for the FFT calculation of the bottom plane currents
fftsizetop=512;

for i=1:size(freq,2)

Jxfftopb(:,:,i)=(fftshift(fft2(Jxb1(:,:,i),fftsizetop,fftsizetop))); %% this is 2D spatial FFT
Jyfftopb(:,:,i)=(fftshift(fft2(Jyb1(:,:,i),fftsizetop,fftsizetop))); 


Jxfftopb2(:,:,i)=(fftshift(fft2(Jxb2(:,:,i),fftsizetop,fftsizetop))); %% this is 2D spatial FFT
Jyfftopb2(:,:,i)=(fftshift(fft2(Jyb2(:,:,i),fftsizetop,fftsizetop))); 


Jmagfftb(:,:,i)=sqrt(abs(Jxfftopb(:,:,i)).^2+abs(Jyfftopb(:,:,i)).^2);
Jmagfftb2(:,:,i)=sqrt(abs(Jxfftopb2(:,:,i)).^2+abs(Jyfftopb2(:,:,i)).^2);


power2d1b(:,:,i)=(Jmagfftb(:,:,i).^2)./max(Jmagfftb(:,:,i).^2,[],'all');
power2d2b(:,:,i)=(Jmagfftb2(:,:,i).^2)./max(Jmagfftb2(:,:,i).^2,[],'all');


end

for i=1:size(Jchangebottom,2)
    Jxfftopb3(:,:,i)=(fftshift(fft2(Jxb3(:,:,i),fftsizetop,fftsizetop))); %% this is 2D spatial FFT
    Jyfftopb3(:,:,i)=(fftshift(fft2(Jyb3(:,:,i),fftsizetop,fftsizetop))); 
    Jmagfftb3(:,:,i)=sqrt(abs(Jxfftopb3(:,:,i)).^2+abs(Jyfftopb3(:,:,i)).^2);
    power2d3b(:,:,i)=(Jmagfftb3(:,:,i).^2)./max(Jmagfftb3(:,:,i).^2,[],'all');
end

% for i=1:size(freq,2)
% EnergynonPc(i)=sum((Jmagfft(:,:,i).^2),'all')./(size(Jmagfft,1).*size(Jmagfft,2));
% EnergyPC(i)=sum((Jmagfft2(:,:,i).^2),'all')./(size(Jmagfft2,1).*size(Jmagfft2,2));
% end
delspacing=abs(nodelist(2,1)-nodelist(1,1));
samplingf1=2*pi*1/delspacing;


fx1=(-(fftsizetop)/2:(fftsizetop)/2-1).*samplingf1./(fftsizetop-1);
fy1=(-(fftsizetop)/2:(fftsizetop)/2-1).*samplingf1./(fftsizetop-1);
[Fx_grid1, Fy_grid1] = meshgrid(fx1, fy1);


%% This portion is for resampling the fft data radially in bins
Radmatrix=sqrt(Fx_grid1.^2+Fy_grid1.^2);

num_bins = 256;  % adjust as needed
r_max = max(Radmatrix(:));
r_edges = linspace(0, r_max, num_bins+1);
% r_centers = (r_edges(1:end-1) + r_edges(2:end)) / 2;

radial_profile = zeros(1, num_bins);
ring_data1 = cell(num_bins, 8); 
ring_data2 = cell(num_bins, 8); 
ring_data3 = cell(num_bins, 8); 

ring_data1b = cell(num_bins, 8); 
ring_data2b = cell(num_bins, 8); 
ring_data3b = cell(num_bins, 8);  
ll=0;

%% This portion samples the fft data radially and then calculates the total power ratio, variance, coefficeint of variation etc
for k=1:13:size(freq,2)
    ll=ll+1;
    temp1=power2d1(:,:,k);
    temp2=power2d2(:,:,k);
    temp3=power2d3(:,:,ll);

    temp1b=power2d1b(:,:,k);
    temp2b=power2d2b(:,:,k);
    temp3b=power2d3b(:,:,ll);
    sumcells(ll)=0;
for i = 1:num_bins
    ring_mask = (Radmatrix >= r_edges(i)) & (Radmatrix < r_edges(i+1));
    
    ring_values1 = temp1(ring_mask);
    ring_values2 = temp2(ring_mask);
    ring_values3 = temp3(ring_mask);

    ring_values1b = temp1b(ring_mask);
    ring_values2b = temp2b(ring_mask);
    ring_values3b = temp3b(ring_mask);

    ring_data1{i,ll}=ring_values1;
    ring_data2{i,ll}=ring_values2;
    ring_data3{i,ll}=ring_values3;

    ring_data1b{i,ll}=ring_values1b;
    ring_data2b{i,ll}=ring_values2b;
    ring_data3b{i,ll}=ring_values3b;

    % sumcells(ll)=sumcells(ll)+size(ring_data1{i,ll},1);
    % scatter(Fx_grid1(ring_mask),Fy_grid1(ring_mask));
    % hold on;
    s1(i,ll)=sum(ring_data1{i,ll});
    s2(i,ll)=sum(ring_data2{i,ll});
    s3_ch5(i,ll)=sum(ring_data3{i,ll});
    s1b(i,ll)=sum(ring_data1b{i,ll});
    s2b(i,ll)=sum(ring_data2b{i,ll});
    s3b_ch5(i,ll)=sum(ring_data3b{i,ll});
end
end
d1=sum(s1,2); %% This is for the total power ratio _top plane
d2=sum(s2,2);
d3_ch5=sum(s3_ch5,2);

d1b=sum(s1b,2);%% This is for the total power ratio _bottom plane
d2b=sum(s2b,2);
d3b_ch5=sum(s3b_ch5,2);

cov1=std(d1(:,1:end),0,2)./mean(d1(:,1:end),2); %% This is for the coefficient of variation_top plane
cov2=std(d2(:,1:end),0,2)./mean(d2(:,1:end),2);
cov3_ch5=std(d3_ch5(:,1:end),0,2)./mean(d3_ch5(:,1:end),2);


cov1b=std(d1b(:,1:end),[],2)./mean(d1b(:,1:end),2); %% This is for the coefficient of variation in the bottom plane
cov2b=std(d2b(:,1:end),[],2)./mean(d2b(:,1:end),2);
cov3b_ch5=std(d3b_ch5(:,1:end),[],2)./mean(d3b_ch5(:,1:end),2);

var1=var(s1,[],2);  %% This is for the variance calculation in top plane
var2=var(s2,[],2);
var3_ch4=var(s3_ch5,[],2);

var1b=var(s1b,[],2);%% This is for the variance calculation in bottom plane
var2b=var(s2b,[],2);
var3b_ch4=var(s3b_ch5,[],2);

%% This part is for plotting the results from the fft analysis


plot(r_edges(2:end),0.5.*(d1+d1b),r_edges(2:end),0.5.*(d2+d2b),r_edges(2:end),0.5.*(d3_ch5+d3b_ch5),r_edges(2:end),0.5.*(d3_ch5+d3b_ch5),r_edges(2:end),0.5.*(d3_ch5+d3b_ch5),r_edges(2:end),0.5.*(d3_ch5+d3b_ch5),r_edges(2:end),0.5.*(d3_ch5+d3b_ch5));
xlim([min(r_edges(2:end)), max(r_edges(2:end))]);
legend('Original','PC1','Design_1','Design_2','Design_3','Design_4','Design_5');
figure;
plot(r_edges(2:end),cc1,r_edges(2:end),cc2,r_edges(2:end),cc3_ch1,r_edges(2:end),cc3_ch2,r_edges(2:end),cc3_ch3,r_edges(2:end),cc3_ch4,r_edges(2:end),cc3_ch5);
xlim([min(r_edges(2:end)), max(r_edges(2:end))]);
legend('Original','PC1','Design_1','Design_2','Design_3','Design_4','Design_5');


imagesc(nodelisthi3t(:,2),nodelisthi3t(:,1),((Jmagfft(9:66,9:66,100))));
xlabel('yaxis');
ylabel('xaxis');
figure;
subplot(3,1,1);
imagesc(nodelisthi3t(:,2).*1e3,nodelisthi3t(:,1).*1e3,((Jmag(15:80-15,15:80-15,51))));
colormap jet;
colorbar;
xlabel('x-axis(mm)',FontSize=14,FontWeight='bold');
ylabel('y-axis(mm)',FontSize=14,FontWeight='bold');
title('|J_{mag,Original}|',FontSize=14,FontWeight='bold')
 set(gca,'YDir','normal');
 rectangle('Position', [-25 -20 50 40], 'EdgeColor', 'b', 'LineWidth', 2); % Red rectangle
subplot(3,1,2);
imagesc(nodelisthi3t(:,2).*1e3,nodelisthi3t(:,1).*1e3,((Jmag2(15:80-15,15:80-15,51))));
colormap jet;
colorbar;
xlabel('x-axis(mm)',FontSize=14,FontWeight='bold');
ylabel('y-axis(mm)',FontSize=14,FontWeight='bold');
title('|J_{mag,Reconstructed}|',FontSize=14,FontWeight='bold');
 set(gca,'YDir','normal');
 rectangle('Position', [-25 -20 50 40], 'EdgeColor', 'b', 'LineWidth', 2); % Red rectangle
 subplot(3,1,3);
imagesc(nodelisthi3t(:,2).*1e3,nodelisthi3t(:,1).*1e3,((Jmag(15:80-15,15:80-15,51)-Jmag2(15:80-15,15:80-15,51))));
colormap jet;
colorbar;
xlabel('x-axis(mm)',Fontsize=14,FontWeight='bold');
ylabel('y-axis(mm)',Fontsize=14,FontWeight='bold');
 set(gca,'YDir','normal');
 title('|J_{mag,Orig}|-|J_{mag,Recon}|',fontsize=14,FontWeight='bold');
 rectangle('Position', [-25 -20 50 40], 'EdgeColor', 'b', 'LineWidth', 2); % Red rectangle





%% This portion is for chromaticity calculation and similarity index calculation
[Gradth2,Gradph2,Gradz]=gradient(ffnorm3,1,1,1);
chrom=chromaticity1(Gradz);
j=0;
for i=1:2:100
    j=j+1;
    % euclideanval(i)=sum(abs((ffnorm2(:,:,i)-ffnorm2(:,:,i+1))).^2,'all');
     [ssimval,ssimmap]=ssim(ffnorm(:,:,i),ffnorm(:,:,i+2),'Radius',0.5);  %% when calculating similarity index for Gradient(Gradients) it is better to scale the values by 1e6. So, the similarity index will bring out the differences.
     % diff(:,:,i)=(Jmag(9:66,9:66,i)./max(Jmag(9:66,9:66,i),[],'all')-Jmag2(9:66,9:66,i)./max(Jmag2(9:66,9:66,i),[],'all'));
      ssval1(j)=ssimval;
      ssmap1(:,:,j)=ssimmap;
end
xvalues=unique(nodelisthi3t(1:end,1));
yvalues=unique(nodelisthi3t(1:end,2));

% This is for making quiver plots
f=51;
figure;quiver(yvalues,xvalues,real(Jcomplex(15:66,15:66,f)),imag(Jcomplex(15:66,15:66,f)),'LineWidth',1.5,'Color','k','AutoScaleFactor',2); %% the x and y are 
hold on;contour(yvalues,xvalues,abs(Jcomplex(15:66,15:66,f)),'LineWidth',1.5);
% set(gca,'YDir','reverse');  %% this shifts the origin to the top left instead of bottom left. 
rectangle('Position', [-0.025 -0.020 0.05 0.04], 'EdgeColor', 'b', 'LineWidth', 2); % Red rectangle
f=1;
figure;
% subplot(1,2,1);
quiver(xvalues, yvalues, real(Jplot2(:,:,f)), imag(Jplot2(:,:,f)), 'r','AutoScale','on','LineWidth',1.5); hold on;
%contour(xvalues,yvalues,abs(Jplot2(:,:,f)));
quiver(xvalues, yvalues, real(Jplot(:,:,f)), imag(Jplot(:,:,f)), 'b','AutoScale','on','LineWidth',1.5); 
% hold on;contour(xvalues,yvalues,Jmag1(:,:,f)-Jmag2(:,:,f),'LineWidth',2);
hold off;
title('Vector Field Comparison');
legend('Original', 'PC1');
axis([-0.1 0.1 -0.1 0.1]);
axis xy;
rectangle('Position', [-0.07 -0.07 0.14 0.14], 'EdgeColor', 'b', 'LineWidth', 2); % Red rectangle
hold on;
rectangle('Position', [-0.017 -0.017 0.034 0.034], 'EdgeColor', 'g', 'LineWidth', 3); % Red rectangle

figure;
contour(fx1,fy1,(Jmagfft2(:,:,5)),3,'r','LineWidth',1);
hold on;
contour(fx1,fy1,(Jmagfft2(:,:,20)),3,'b','LineWidth',1);
hold on;
contour(fx1,fy1,(Jmagfft2(:,:,60)),3,'g','LineWidth',1);
hold on;
contour(fx1,fy1,(Jmagfft2(:,:,90)),3,'k','LineWidth',1);
axis xy;
rectangle('Position', [-0.07 -0.07 0.14 0.14], 'EdgeColor', 'b', 'LineWidth', 2); % Red rectangle
hold on;
rectangle('Position', [-0.017 -0.017 0.034 0.034], 'EdgeColor', 'g', 'LineWidth', 3); % Red rectangle
title('|Jpc1|-|Joriginal|');
colormap jet;


%% This is used for importing the farfield data for any new designs and immediately evaluating their simialrity index
allfile3=dir('F:\rri files\APSERA\surface current data\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_140mm\Final_codes_with_sphericalsampling\patch_original\Farfield\*[1]*.txt');
for i=1:size(allfile3,1)
    basefile3=allfile3(i).name;
    filena3=fullfile('F:\rri files\APSERA\surface current data\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_140mm\Final_codes_with_sphericalsampling\patch_original\Farfield',basefile3);
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
for k=1:1:size(farff4,3)
    for j=1:size(farff4,2)
            maxval=max(farff4(:,j,k),[],'all');
           ffnorm4(:,j,k)=farff4(:,j,k)./maxval;
    end

end

[Gx,Gy,Gz1]=gradient(ffnorm4(:,:,:),1,1,1);
chrom=chromaticity1(Gz1(:,:,1:end));




%% This portion is for creating a movie from the currents

% Example input: J_all (M x N x 100), complex current matrices
% Customize the filename and frame delay
filename = 'current_animation2.gif';
delayTime = 0.5;  % seconds between frames
diff1=diff(9:66,9:66,:);
M = size(diff1,1);
N = size(diff1,2);

% Normalize across all frames for consistent color scaling
max_val = max((diff1(:)));
min_val = min((diff1(:)));

for k = 1:size(diff1, 3)
    % Extract current frame
    J = diff1(:,:,k);

    % You can use magnitude or phase; here we use magnitude
    imagesc(nodelisthi3t(:,1),nodelisthi3t(:,2),J, [min_val, max_val]);
    xlabel('x-axis');
    ylabel('y-axis');

    % axis image off;
    colormap(jet);
    colorbar;
    title(['Frequency ', num2str(freq(k))]);

    % Capture frame as image
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);

    % Write to GIF
    if k == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', 5, 'DelayTime', delayTime);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end



