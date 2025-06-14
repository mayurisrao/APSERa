clear all;

%% load the appropriate currents file at multiple frequencies.
load ('D:\Keerthi\PMA_tapered_circular patch\current_patch.mat');
% Jmin=min(abs(recon_J),[],1);
% Jmax=max(abs(recon_J),[],1);
% Jt=((re-Jmin)./(Jmax-Jmin));
   Jt=J; 
 nodelist=nodelist;
%% load the farfield data
 load('D:\Keerthi\Final_codes_with_sphericalsampling\farield_test1.mat');
 
freq=(2.5:0.015:4).*1e9;  %% this is for 101 points, this will get changed if number of points are different
%% calculation of similarity index
for i=1:size(freq,2)-1
    % in this line we are copmaring the beam at 3GHz with the rest of the
    % beams  
     % vari_(i)=sum(abs((directivitynorm(:,:,i)-directivitynorm(:,:,12))).^2,'all');
     % if(vari_(i)==0)
     %     vari_(i)=0.01;
     % end
     % vari(i)=1./(vari_(i));
     vari2_(i)=sum(abs((directivitynorm(:,:,i)-directivitynorm(:,:,49))).^2,'all');
     % if(vari2_(i)==0)
     %     vari2_(i)=0.01;
     % end
     if(i>1 && vari2_(i)==0)
        vari2(i)=1/((vari2_(i-1))+vari2_(i));
     else
        vari2(i)=1/vari2_(i);
     end
     
     vari1(i)=sum(abs((directivitynorm(:,:,i)-directivitynorm(:,:,i+1))).^2,'all');
     varwei1(i)=vari2(i)./vari1(i);
     [ssimval,ssimmap]=ssim(directivitynorm(:,:,i).*1e4,directivitynorm(:,:,i+1).*1e4,'Radius',1);   %% the directivity values are scaled by 1e3 because the similarity index recognizes the differences more prominently
     simval(i)=ssimval;
     ssmap(:,:,i)=ssimmap;
end
varwei1(end+1)=varwei1(end);

Jxtop=Jt(1:size(nodelist,1),:);
Jytop=Jt(size(nodelist,1)+1:end,:);
% Jxtop=Jt(1:4096,:);
% Jytop=Jt(4097:end,:);
Jsxtshaped=reshape(Jxtop,sqrt(size(nodelist,1)),sqrt(size(nodelist,1)),51);  %% incase of different number of frequency points, change the last index
Jsytshaped=reshape(Jytop,sqrt(size(nodelist,1)),sqrt(size(nodelist,1)),51);
% Jsxtshaped=reshape(Jxtop,70,70,101);  %% incase of different number of frequency points, change the last index
% Jsytshaped=reshape(Jytop,70,70,101);


% load similarityindex.mat;  %% this has the similarity index 
%% This calculates the weights for the frequencies from the similarity index
% varwei1(1:100)=1+(((simval-min(simval))/(max(simval)-min(simval))).*1000);
%   varwei1(1:30)=0.00001;
%   %varwei1(1:10)=0.01;
% 
%   varwei1(80:101)=0.00001;
% %% inorder to reduce the impact of the the higher edge frequencies, we are assigning the last few weights to 0.01.
 % varwei1(1:51)=1;
% performing non centered weighted pca

[Jxcoefftop,Jxscoretop,explainedxtop,mux]=eigenowncode(Jxtop,varwei1);

[Jycoefftop,Jyscoretop,explainedytop,muy]=eigenowncode(Jytop,varwei1);

%% reconstructing the currents with the PCA
% Jxcoefftop=10.^(Jxcoefftop);
% Jycoefftop=10.^(Jycoefftop);
% Jxscoretop=10.^(Jxscoretop);
% Jyscoretop=10.^(Jyscoretop);
reconxtop=[1:1];  %% specify the PCA component indices here for the Jx component
reconytop=[1:1];  %% specify the PCA component indices here for the Jy component
recon_jxtop=Jxscoretop(:,reconxtop)*Jxcoefftop(:,reconxtop)';
recon_jxtop=((recon_jxtop+mux))./varwei1;

recon_jytop=Jyscoretop(:,reconytop)*Jycoefftop(:,reconytop)';
recon_jytop=((recon_jytop+muy))./varwei1;
  
recon_J=[recon_jxtop;recon_jytop];

% recon_Jbottom=recon_J;
%% this saves the reconstructed current in a file
      save('D:\Keerthi\PMA_tapered_circular patch\J_top_reconv1.mat','recon_J','nodelist','Jxcoefftop','Jycoefftop','Jxscoretop','Jyscoretop','varwei1');

% reconJx=reshape(recon_J(1:4096,:),64,64,101);
% reconJy=reshape(recon_J(4097:end,:),64,64,101);
% reconJJ=sqrt(abs(reconJx).^2+abs(reconJy).^2);