clearvars;
close all;
%% read all the files from the E-field folder
allfile=dir('D:\Keerthi\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_surfacecurrents - slots22\Export\3d\*[1].txt*');
for i=1:size(allfile,1)
    basefile1=allfile(i).name;
    filena1=fullfile('D:\Keerthi\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_surfacecurrents - slots22\Export\3d',basefile1);
    pattern = 'f=(\d+\.\d+|\d+)'; %% this is to identify the frequency values automatically from the file names
    tokens = regexp(filena1, pattern, 'tokens');
    if ~isempty(tokens)
    freq(i) = str2double(tokens{1}{1});
    disp('Extracted frequency:');
    disp(freq(i)); 
    end 
    c1(:,:,i)=importdata(filena1).data;
end


freq=sort(freq,'ascend');
 %% determine the x and y values
xloc=(c1(:,1,1)-70).*1e-3; %% this is for having the x and y centered at 0,0 %% this should be changed according to the plane
yloc=(c1(:,2,1)-70).*1e-3;%% this is for having the x and y centered at 0,0 %% this should be changed according to the plane dimensions
zloc=(c1(:,3,1)).*1e-3;
zvalues=(unique(c1(:,3,1)) ).*1e-3;
yvalues=((unique(c1(:,2,1)-70))).*1e-3;
xvalues=((unique(c1(:,1,1))-70)).*1e-3;

%%this was the set of points used for sampling in CST
theta=-90:2:90;
phi=0:2:180;
r=150;  %% units in mm.( this is the format CST takes)
k=0;
for i=1:size(phi,2)
    for j=1:size(theta,2)
        k=k+1;
        sample(k,1)=r*sind(theta(j))*cosd(phi(i));
        sample(k,2)=r*sind(theta(j))*sind(phi(i));
        sample(k,3)=r*cosd(theta(j));
    end
end
sample(:,1)=(sample(:,1)+70);
sample(:,2)=(sample(:,2)+70);
sample(:,3)=(sample(:,3)+2.1);

[p,t]=meshgrid(phi,theta);
ptop=p(:);
ttop=t(:);

efieldxcomplextop=zeros(size(c1,1),size(freq,2));
efieldycomplextop=zeros(size(c1,1),size(freq,2));
efieldzcomplextop=zeros(size(c1,1),size(freq,2));
efieldvectortop=zeros(size(c1,1),3,size(freq,2));

% efieldxtopshaped=zeros(size(xval,2),size(yval,2),size(freq,2));
% efieldytopshaped=zeros(size(xval,2),size(yval,2),size(freq,2));
% efieldztopshaped=zeros(size(xval,2),size(yval,2),size(freq,2));

%% extracting the E-field top
for i=1:size(allfile,1)
    efieldxcomplextop(:,i)=c1(:,4,i)+1i*c1(:,5,i);
    efieldycomplextop(:,i)=c1(:,6,i)+1i*c1(:,7,i);
    efieldzcomplextop(:,i)=c1(:,8,i)+1i*c1(:,9,i);
    efieldvectortop(:,:,i)=[efieldxcomplextop(:,i) efieldycomplextop(:,i) efieldzcomplextop(:,i)];
%     efieldxtopshaped(:,:,i)=transpose(reshape(efieldxcomplextop(:,i),size(xval,2),size(yval,2)));
%     efieldytopshaped(:,:,i)=transpose(reshape(efieldycomplextop(:,i),size(xval,2),size(yval,2)));
%     efieldztopshaped(:,:,i)=transpose(reshape(efieldzcomplextop(:,i),size(xval,2),size(yval,2)));

end


%% creation of nodes in a rectangular plane for the surface curretnts
xrange=0.205;   %% this range should ideally be kept 50% more than the antenna size
yrange=0.205;     %% this range should ideally be kept 50% more than the antenna size
Nx=75;Ny=75;    %% this defines the number of nodes in the eq.current plane. eg. For Nx=71, there will be Nx-1 (70) nodes.
hx=xrange/(Nx-1); 
hy=yrange/(Ny-1); 
cellarea=hx*hy;
% hx=3.5;hy=3.5; %% this can be changed according to the required sampling interval
totalnodes=(Nx-1)*(Ny-1); 
nodelist=zeros(totalnodes,3); 
k=1;
zs=0.000 ;  %% this defines the z axis location of the equivalent current (this should be zero as we have assumed the equivalent surface current plane being at origin)

%% calculation of nodes in source current plane
for ll=1:size(zs,2)
for i=1:Ny-1
        for j=1:Nx-1
            xi=(-xrange/2)+(hx/2)+((j-1)*hx);
            yi=(-yrange/2)+(hy/2)+((i-1)*hy);
%             l=xmin+(j-1)*hx;
%             m=ymin+(i-1)*hy;
            nodelist(k,:)=[xi yi zs];
            k=k+1;
        end
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
end



%% use the below lines to write the sample points to a text file in a format compatible to CST. 
% 
%  fileID = fopen('sample_dipoletop.txt', 'w');  % 'w' for writing (creates or overwrites the file)
% 
% % Write the matrix to the file, separating columns by tabs
% for i = 1:size(sample, 1)  % Loop through each row
%     fprintf(fileID, '%f\t%f\t%f\n', sample(i, :));  % Print each row with tabs separating columns
% end
% 
% % Close the file
% fclose(fileID);

% Close the file
% fclose(fileID);

%% converting the cartesian Efield to spherical E field     
Ethetatopf=zeros(size(theta,2),size(phi,2),size(freq,2));
Ephitopf=zeros(size(theta,2),size(phi,2),size(freq,2));
Ettopf=zeros(size(xloc,1),size(freq,2));
Eptopf=zeros(size(xloc,1),size(freq,2));
for k=1:size(freq,2)
    ll=0;
    for i=1:size(phi,2)
        for j=1:size(theta,2)
            ll=ll+1;
            Ethetatopf(j,i,k)=[cosd(theta(j)).*cosd(phi(i)) cosd(theta(j)).*sind(phi(i)) -sind(theta(j))]*[efieldxcomplextop(ll,k);efieldycomplextop(ll,k) ;efieldzcomplextop(ll,k)];
            Ephitopf(j,i,k)=[-sind(phi(i)) cosd(phi(i)) 0]*[efieldxcomplextop(ll,k);efieldycomplextop(ll,k) ;efieldzcomplextop(ll,k)];
        end
    end
    Ettopf(:,k)=reshape(Ethetatopf(:,:,k),size(xloc,1),1);
    Eptopf(:,k)=reshape(Ephitopf(:,:,k),size(xloc,1),1);

end

w=2.*pi.*freq.*10^9;
mu=4*pi*(1e-7);
epsi=8.854*1e-12;
eeta=377;
lamda0=3e8./freq/1e9;
beta0=2*pi./lamda0;
Q11t=zeros(size(xloc,1),size(nodelist,1));
Q12t=zeros(size(xloc,1),size(nodelist,1));
Q21t=zeros(size(xloc,1),size(nodelist,1));
Q22t=zeros(size(xloc,1),size(nodelist,1));

fieldarraytop=[xloc yloc zloc];
diff_matrix = reshape(fieldarraytop, [size(fieldarraytop, 1), 1, 3]) - reshape(nodelist, [1, size(nodelist, 1), 3]);
diff_matrix2=(diff_matrix.^2);
dist = (sqrt(sum(diff_matrix2, 3))); % calculates Pairwise distances
dist2=dist.^2;
dist3=dist2.*dist;
dist4=dist3.*dist;
dist5=dist4.*dist;
n=size(nodelist,1);

for k=1:size(freq,2)  % this loop defines the frequency points. This can be changed according to the range of calculation (i.e either one frequency or multiple frequency point calcs)
    B=[Ettopf(:,k);Eptopf(:,k)];  %% this calculates the RHS of the MOM eqn
    const1=(1j*w(k)*mu./(dist3))+(3*eeta./(dist4))+(3./(1j.*w(k).*epsi.*(dist5)));
    
    const2=(1j*w(k)*mu./dist)+(eeta./(dist2))+(3./(1j*w(k)*epsi.*(dist3)));
    const3=exp(-1j.*beta0(k).*dist);
    const4=(1j*w(k)*mu./(dist3))+(3*mu./(dist4))+(3./(1j.*w(k).*epsi.*(dist)));
    const5=(1j*w(k)*mu./(dist3))+(3*mu./(dist4))+(3./(1j.*w(k).*epsi.*(dist5)));
    const6=(1j*w(k)*mu./dist)+(mu./(dist2))+(1./(1j*w(k)*epsi.*(dist3)));
    for i=1:size(xloc,1)
        for j=1:size(nodelist,1)

            constant1=const1(i,j);
            constant2=const2(i,j);
            constant3=const3(i,j);
            constant4=const4(i,j);
            constant5=const5(i,j);
            constant6=const6(i,j); 
            dx=diff_matrix(i,j,1);
            dy=diff_matrix(i,j,2);
            dz=diff_matrix(i,j,3);
            dx2=diff_matrix2(i,j,1);
            dy2=diff_matrix2(i,j,2);
            dxdy=dx*dy;
            convercon1=cosd(ttop(i)).*cosd(ptop(i));
            convercon2=cosd(ttop(i)).*sind(ptop(i));
            convercon3=sind(ttop(i));
            Q11t(i,j)=((convercon1*constant3*((dx2*constant1)-constant2))+(convercon2*constant3*(dxdy)*constant1)-(convercon3*constant3*(dx)*(dz)*constant1))*cellarea/4/pi;
            Q12t(i,j)=((convercon1*constant3*(dxdy)*constant1)+(convercon2*constant3*((dy2*constant1)-constant2))-(convercon3*constant3*(dy)*(dz)*constant1))*cellarea/4/pi;
            Q21t(i,j)=((-cosd(ptop(i))*constant3*(dxdy)*constant1)-(sind(ptop(i))*constant3*((dx2*constant1)-constant2)))*cellarea/4/pi;
            Q22t(i,j)=((cosd(ptop(i))*constant3*((dy2*constant5)-constant6))-(sind(ptop(i))*constant3*(dxdy)*constant1))*cellarea/4/pi;
            
        end
    end

   H=[Q11t Q12t;Q21t Q22t];  %% this is the matrix to be inverted to solve for the currents
   [J(:,k),flag(k),res(k)]=lsqr(H,B,1e-5,3000);  %% this will also store the flag and the residual error after solving.
   disp(k);
   disp(rms(H*J(:,k)-B));
   if(flag(k)==0)
       disp('converged'); 
   else
       disp('not converged');
   end
    
end 
Jhi3=J;
% Jtemp=J(:,50);

nodelisthi3=nodelist;
save('D:\Keerthi\DIPOLE_2GHZ_dualpol_TEST_broadband_WITH slot_9_with sma_finitegroundplane3_surfacecurrents - slots22 \Export\3d\current.mat','Jhi3','nodelisthi3'); 

  % Jx1_=reshape(Jhi3(1:5476,:),74,74,8);
  % Jy1_=reshape(Jhi3(5477:end,:),74,74,8);
%   nodelistshapedx=reshape(nodelist(:,1),64,64);
%   nodelistshapedy=reshape(nodelist(:,2),64,64); 
%   nodelisthi3shapedx=reshape(nodelisthi3(:,1),64,64);
%   nodelisthi3shapedy=reshape(nodelisthi3(:,2),64,64);
% 
% for i=1:20
%     Jx(:,:,i)=interp2(nodelisthi3shapedx',nodelisthi3shapedy',Jx1_(:,:,i),nodelistshapedx',nodelistshapedy',"nearest");
%     Jy(:,:,i)=interp2(nodelisthi3shapedx',nodelisthi3shapedy',Jy1_(:,:,i),nodelistshapedx',nodelistshapedy','nearest');
%     Jx_(:,i)=reshape(Jx(:,:,i),4096,1);
%     Jy_(:,i)=reshape(Jy(:,:,i),4096,1);
%     Jt_(:,i)=[Jx_(:,i);Jy_(:,i)];
% end
% 
% 
% load('D:\Keerthi\Final_codes_with_sphericalsampling\current_total2.mat');
% 
% 

% 
% 
% 
