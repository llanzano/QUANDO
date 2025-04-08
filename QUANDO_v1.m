%%   User-fiendly code based on the following article:
%$£   
%$£  "Location of oncogene-induced DNA damage sites revealed by quantitative analysis of a DNA counterstain"
%$£  by Greta Paterno' et al , European Biophysics Journal 2025
%$£
%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£ 
%$£                                                                     %$£
%$£                              Luca Lanzanò                           %$£
%$£       University of Catania - Department of Physics and Astronomy   %$£
%$£           & Istituto Italiano di Tecnologia - Nanoscopy             %$£
%$£                      User-Friendly Version (........)               %$£
%$£                                                                     %$£
%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£%$£


function [  ]=QUANDO_v1(filename );

% OUTPUT text FILE contain: 
% First 3 columns:
% cell analyzed, DNA density, Colocalization with Heterochromatin
% The last 3 columns: 
% file analyzed, frame analyzed, cell # in the count mask
%
% goal of the script:
%
% for 2 channels
% -normalize intensity of DNA marker in each cell (1 max value)
% -estimate the average value in the nucleus and in the loaded ROI mask
% corresponding to DNA damage foci

% load: 
% count mask nuclei 
% ROI mask 
% intensity image (specify number of channels and DNA channel)
%(NOTE: AVOID SALT and PEPPER NOISE. If Present, denoise image with
% imagej using 'process'--> 'noise' --> 'despeckle')

% open stacks

% defaults values TO CHECK
cellfirst=1;
cellmax=5000;

DNAnormTHR=0.3;
nch=2;
chDNA=2; % channel for DNA (DAPI)
tag=['Results'];
ThrImg1=6000;
ThrImg2=1000;
maxlag=14;

% menu for channels
prompt = {'Number of channels', 'DNA channel', 'Output filename', 'ICCS Threshold ch1', 'ICCS Threshold DNA channel', 'ICCS max lag'  }; 
dlg_title = 'Press OK to confirm'  ; 
num_lines = [1 50];
def = { num2str(nch), num2str(chDNA), tag, num2str(ThrImg1), num2str(ThrImg2), num2str(maxlag)  };
answer = inputdlg(prompt,dlg_title,num_lines,def);
nch=str2num(answer{1});
chDNA=str2num(answer{2});
tag=answer{3};
ThrImg1=str2num(answer{4});
ThrImg2=str2num(answer{5});
maxlag=str2num(answer{6});

ch1=1;
if chDNA==1 && nch==2
ch1=2 ;
end

mask=0;  % mask to use as ROI 
maskedu=1;
doublech=0; % flag for 1ch or 2ch analysis
maskinput=1;  % load n number of masks
%insert loop for multiple files
% if isfile('list.mat')
%     load('list.mat')
% else
[List]=ICCS_SelectFiles(maskinput) ; % description: load Count mask nuclei, (copy), intensity file, ROI Mask
save('list.mat', 'List' );
% end
cell=0;
% spot=0;
% spot2=0;
% dist=0;
% ResultSpot1=double.empty(0,9);
% ResultSpot2=double.empty(0,8);
ResultCell=double.empty(0,16);
ResultCellSel=double.empty(0,6);
% ResultCellG1=double.empty(0,11);
% ResultCellEarly=double.empty(0,11);
% ResultCellLate=double.empty(0,11);
% ResultDist=double.empty(0,7);

for ifile=1:size(List,1)
%open files ... 
% count mask nuclei
filenamefull1=List{ifile,1};
MCstack=simple_ICCS_readfiletif(filenamefull1);  % read stack of count mask nuclei
filenamefulle = List{ifile,2};
MCstack1=simple_ICCS_readfiletif(filenamefulle); % read stack of count mask spots (nuclei in this file)
if size(List,2)>2 &&~isempty(List{ifile,3})
filenamefulle2 = List{ifile,3};
Astack=simple_ICCS_readfiletif(filenamefulle2);  % order :  EdU, ch2, Dapi
end
if maskedu>0 
filenamefulle4 = List{ifile,4};
MaskEduStack=simple_ICCS_readfiletif(filenamefulle4);  % load Mask for analysis of EdU pixels
end
for istack=1:size(MCstack,3)
MC=MCstack(:,:,istack);
MinN=min(MC(MC>0),[], 'all');
MaxN=max(MC,[], 'all');
% count mask CH1
MC1=MCstack1(:,:,istack);
% intensity file
if size(List,2)>2 &&~isempty(List{ifile,3})
A=Astack(:,:,(istack-1)*nch+1:(istack-1)*nch+nch);
end
if maskedu>0 
MaskEdu=MaskEduStack(:,:,istack);
else
MaskEdu=MC;  % 
end

if mask>0 && ~isempty(List{ifile,3+mask})  %mask>0?  &&   OK per ora
filenamefulle4 = List{ifile,3+mask};
MaskROI=simple_ICCS_readfiletif(filenamefulle4);  % load Mask for ROI of analysis
else
MaskROI=MC;  % 
end

X=size(MC,1);
Y=size(MC,2);
MinCell=max([cellfirst, MinN]);
MaxCell=min([cellmax, MaxN]);
Imgall=zeros(size(MCstack(:,:,1),1),size(MCstack(:,:,1),2));
for cellidx=MinCell:MaxCell
cell=cell+1;
% mask of single cell (value 1)    
MaskRaw=zeros(X,Y);
MaskRaw(MC==cellidx)=1;
% MaskRaw(MaskROI==0)=0;  % apply mask defining ROI, if loaded
%mask of spots in ch1 (values from k to k+nspots)
Mask1=MC1;
Mask1(MaskRaw==0)=0;

% operate on DNA image: must normalize for each cell
Imgdapi=double(A(:,:,chDNA));
Imgdapi(MaskRaw==0)=0;
% calculate intensities
Imgch1=double(A(:,:,ch1));
Imgch1(MaskRaw==0)=0;
[ Mraw1, Mraw2 ]=ICCS_for_QUANDO( Imgch1, Imgdapi, MaskRaw, ThrImg1, ThrImg2, maxlag ) ; 

% Imgdapi=imtranslate(Imgdapi,[-120/45,-40/45]);     %translation
%NORMALIZATION
Imgdapinorm=Imgdapi./max(Imgdapi,[],'all');  % normalized dna image for each cell
Idapinorm=mean(Imgdapinorm(MaskRaw>0));
IdapinormStDev=std(Imgdapinorm(MaskRaw>0));
Imgdapimask=Imgdapinorm; 
Imgdapimask(Imgdapinorm>=DNAnormTHR)=1; % mask for high DAPI signal 
Imgdapimask(Imgdapinorm<DNAnormTHR)=0;

Idapinormmask=0;
HCI=0;
I1mask=0;
Idapimask=0;
% Idapistdmask=0;
CVmask=0;
% Ispot=0;

Npxdapi = numel(find( MaskRaw>0 )) ;
Idapi=mean(Imgdapi(MaskRaw>0));

Npxedu = numel(find(MaskEdu>0 & MaskRaw>0 )) ; % pixels in loaded mask (e.g. edu mask)
if Npxedu>0
Idapinormmask=mean(Imgdapinorm(MaskEdu>0 & MaskRaw>0));  % Idapinormmask = av value of normalized DNA image in mask pixels !!!
NpxeduHC = numel(find(MaskEdu>0 & MaskRaw>0 & Imgdapimask>0 )) ;  % pixels with high DAPI signal 
HCI=NpxeduHC/Npxedu; % HC-index = percentage of pixels with high DAPI signal

I1mask=mean(Imgch1( (MaskEdu>0 & MaskRaw>0)  ));
Idapimask=mean(Imgdapi( (MaskEdu>0 & MaskRaw>0)  ));
Idapistdmask= std(Imgdapi( (MaskEdu>0 & MaskRaw>0)  ));
CVmask =  Idapistdmask / Idapimask ;

end
% Idapinormmask=NpxeduHC/Npxedu;
% if N1>0        
% Ispot=mean(Imgspot(Mask1>0));  
% end
Imgall=Imgall+Imgdapinorm;
ResultsCellTemp=cat(2, cell , Idapinorm, IdapinormStDev, Idapinormmask, HCI, Npxdapi, Idapi, Npxedu,  I1mask,  Idapimask ,  CVmask , Mraw1, Mraw2, ifile, istack, double(cellidx)   )  ;
%extract for ech cell: 
ResultsCellTempSel=cat(2, cell , Idapinormmask, Mraw1, ifile, istack, double(cellidx)   )  ;

% ResultsCellTemp=cat(2, cell , N1,  SizeAv1, Ispot, Npxdapi, Idapi, Npxedu, Iedu, ifile, istack, double(cellidx)   )  ;
% ResultsCellTemp=cat(2, cell , N1,  SizeAv1, N2,  SizeAv2, Ncoloc , ifile, double(cellidx)   )  ;
ResultCell=cat(1,ResultCell,ResultsCellTemp) ; % must initialize
ResultCellSel=cat(1,ResultCellSel,ResultsCellTempSel) ; % must initialize
 
close all
% end
end
% figure
% imagesc(Imgall)
% NormDNA(istack)=Imgall;
end


end

% filenameout=[tag ];   % cell #,  ... 
% dlmwrite([filenameout,'.txt'],ResultCell,'delimiter',';','precision',8);
filenameout=[tag ];   % cell #,  ... 
dlmwrite([filenameout,'.txt'],ResultCellSel,'delimiter',';','precision',8);

end


%required funtions

function A=simple_ICCS_readfiletif(fname)
info = imfinfo(fname);
nslice = numel(info);
A=imread(fname, 1);  % read tif stack data
for k = 2:nslice
    B = imread(fname, k);
    A=cat(3,A,B);
end

end

function [List]=ICCS_SelectFiles(mask) ;

%open files ... 
filename='none';
c=1;
while filename~= 0
[filename,pathname, filterindex] = uigetfile({'*.tif'},['Select Count Mask Nuclei ',num2str(c)]);
if filename~= 0
filenamefull = [pathname, filename];  
List{c,1}=filenamefull;

% [filenamemask,pathnamemask, filterindex] = uigetfile({'*.tif'},['Select Count Mask Nuclei ',num2str(c)]);
% filenamefull3 = [pathnamemask, filenamemask]; 
% List{c,2}=filenamefull3;
List{c,2}=filenamefull;

[filenamemaske,pathnamemaske, filterindex] = uigetfile({'*.tif'},['Select Intensity file (Despeckle) ',num2str(c)]);
if filenamemaske ~= 0
filenamefulle = [pathnamemaske, filenamemaske];   
List{c,3}=filenamefulle ;
end
if mask>0
    for imask=1:mask
    [filenamemaskr,pathnamemaske, filterindex] = uigetfile({'*.tif'},['Select ROI Binary Mask',num2str(imask),' of file ',num2str(c)]);
        if filenamemaskr ~= 0
        filenamefullr = [pathnamemaske, filenamemaskr];   
        List{c,3+imask}=filenamefullr ;
        end
    end
end

c=c+1;
end

end
end


function [ Mraw1, Mraw2 ]=ICCS_for_QUANDO( xch1, xch2, MaskRaw, ThrImg1, ThrImg2, maxlag )

% defaults values
save_each_cell=0; % 1: save the analysis for each cell. 0: only AllResults-ICCS file
ch2=2;  % 2nd channel for analysis! write: 2 or 3 or '2+3'
menuflag=0;
FWHM1nm=200;
FWHM2nm=200;
px=0.045; %px size in um
% ThrImg1=150;    % SET Threshold !!! (non necessario con EdU  mask)
% ThrImg2=500;   % SET Threshold !!!
Sat1=0;    % saturate intensity above this value (if >0)
Sat2=0; % saturate intensity above this value (if >0)
% lag00=2; %first spatial lag to fit
lag0=[0 1 1];
% maxlag=14;
nongaussian=0;  % 0 gaussian, 1 use new formula with cosh
M=10;  % points in x for local analysis (1 for global only)- do not use in this script
Thrmask=0;
MinPxPos=0; % min pixels in 2nd mask to detect positive cell (0 for analyzing ALL cells)
cell0=0; % cell to visualize for checking fit (0 for none)
cellfirst=1;
cellmax=10000000;  
AllResults=double.empty(21,0);

%insert loop for multiple files
% [List]=ICCS_SelectFiles() ; %
% save('list-iccs.mat', 'List' );
   
%set analysis area (if no mask is loaded)
Thr1=-1;
Thr2=-1;
sm=0.0;
Extra=0;
figcheck=1;
% maxlag=ceil(X/4);
figure

% B=simple_ICCS_Threshold(x1,Thr1,Thr2);
xch1t=xch1;
xch2t=xch2;
xch1t(xch1<ThrImg1)=0;  % this is based on numerical common threshold
% xch1t(Mask2ndRaw==0)=0; % this is based on loaded EdU Mask
xch2t(xch2<ThrImg2)=0;
x1=cat(3,xch1t,xch2t);


Mask=simpleICCS_Threshold(MaskRaw,Thrmask);
[x1pad,Np,Mask,B,Iav]=simple_PadForICS_fromMask(x1,Extra,Mask);

b1=x1pad(:,:,1);
[m,n]=size(b1);
b2=x1pad(:,:,2);
Iav(1)=mean(b1(Mask>0));
Iav(2)=mean(b2(Mask>0));

xch1tm=xch1t;
xch2tm=xch2t;

% subplot(2,2,1)
% imagesc(b1)
% subplot(2,2,2)
% imagesc(b2)
% subplot(2,2,3)
% imagesc(Mask)

%ICS global analysis with padding G(0) correction
ACF1=simple_ICCS_CCFmean(b1,b1)*m*m*n*n./(m*n-Np(1))./(m*n-Np(1));
ACF2=simple_ICCS_CCFmean(b2,b2)*m*m*n*n./(m*n-Np(2))./(m*n-Np(2));
AreaRatio1=( m*n-Np(1) ) / ( m*n ) ;
AreaRatio2=( m*n-Np(2) ) / ( m*n ) ;
CCF12=simple_ICCS_CCFmean(b1,b2)*m*m*n*n./(m*n-Np(1))./(m*n-Np(2));  %Np(1) and 2 are equal
CCF21=simple_ICCS_CCFmean(b2,b1)*m*m*n*n./(m*n-Np(1))./(m*n-Np(2));  %Np(1) and 2 are equal
CCF=0.5*(CCF12+CCF21);

maxlag=min(maxlag,length(ACF1));
% subplot(2,2,4)
% plot(ACF1(1:maxlag),'DisplayName','ACF1');hold all;plot(ACF2(1:maxlag),'DisplayName','ACF2');plot(CCF(1:maxlag),'DisplayName','CCF');hold off;

%set gamma factor
gamma=0.5;
%set parameters for fit
% lag00=1; %first spatial lag to fit
lagfit0=maxlag; %points to fit
tol=50;
% const=0.19; 
% const=simple_ICCS_Get_const(FWHM1nm/(px*1000),FWHM2nm/(px*1000)); 
% lag0=[lag00 lag00 lag00];
lagfit=[lagfit0 lagfit0 lagfit0];
% nongaussian=0;  % 0 gaussian, 1 use new formula with cosh

    if length(lag0)<3
        if length(lag0)<2
        lag0(2)=lag0(1);
        lag0(3)=lag0(1);
        else
            lag0(3)=lag0(2);
        end
    end
    if length(lagfit)<3
        if length(lagfit)<2
        lagfit(2)=lagfit(1);
        lagfit(3)=lagfit(1);
        else
            lagfit(3)=lagfit(2);
        end
    end
AreaPSFRatio1=m*n/(pi*2*((FWHM1nm/1000/px)*0.25*sqrt(2/log(2)))^2);
AreaPSFRatio2=m*n/(pi*2*((FWHM2nm/1000/px)*0.25*sqrt(2/log(2)))^2);

subplot(2,2,3)
[param1, fit1, chisq1]=simple_ICCS_Fit_ICS_1Dsingle(lag0(2):lag0(2)+lagfit(2)-1,ACF1(lag0(2)+1:lag0(2)+lagfit(2)),4,1,'auto 1');
subplot(2,2,4)
[param2, fit2, chisq2]=simple_ICCS_Fit_ICS_1Dsingle(lag0(3):lag0(3)+lagfit(3)-1,ACF2(lag0(3)+1:lag0(3)+lagfit(3)),4,1,'auto 2');
w1=param1(3);
w2=param2(3);
Size1=w1*px*1.18;
Size2=w2*px*1.18;
wrefg=sqrt(w1*w2);
wrefdelta=sqrt(0.5*( w1^2+w2^2 ));

%other parameters
Br1=param1(2)*AreaRatio1*Iav(1);
N1=(gamma/param1(2))*AreaPSFRatio1;
Br2=param2(2)*AreaRatio2*Iav(2);
N2=(gamma/param2(2))*AreaPSFRatio2;


Mraw1 = CCF(1) *( 1/param2(2)) ; 
Mraw2 = CCF(1) *( 1/param1(2)) ; 
Mraw = 0.5*( Mraw1 + Mraw2 ) ; 


close all

end

%required funtions ICCS

function [y,Np, varargout]=simple_PadForICS_fromMask(x1,Extra,Mask)
[m,n,p]=size(x1);
Mask=double(Mask);
%% padding
y=zeros(m+2*Extra,n+2*Extra,p);
for k=1:p
    y(Extra+1:Extra+m,Extra+1:Extra+n,k)=x1(:,:,k);
end
%% adding average on zeros
for k=1:p
    x=y(:,:,k);
    MeanInt(k)=mean(x(Mask>0));
    c=0;
    for i=1:m+2*Extra
        for j=1:n+2*Extra
            if Mask(i,j)==0 || isnan(Mask(i,j))
            y(i,j,k)=MeanInt(k) ;
            c=c+1;
            end
        end
    end
Np(k)=c;

if nargout > 2
varargout{1} = Mask;
end

if nargout > 3
    for k=1:p
      Aroi=y(:,:,k).*Mask;
      A=simpleICCS_smooth_simple(Aroi,0.2,1);
      B(k)=median(A(A>0));
    end    
varargout{2} = B;
end

if nargout > 4   
varargout{3} = MeanInt;
end



    
end

end


function  Output=simple_ICCS_CCFmean(x1,x2)

NumberOfAngles=180;
[X,Y]=size(x1);
%ACF=conv2(x1,x2,'same');
F1=fft2(x1);
F2=fft2(x2);
ACF= F1.*conj(F2);
G=((sum(sum(x1)))*(sum(sum(x2)))/X/Y);
ACF= ifft2(ACF);
ACF= fftshift(ACF)./G-1;

[R, C]=size(ACF);
if iseven(R)
r0=R/2+1;
else
r0=(R+1)/2;
end
if iseven(C)
c0=C/2+1;
% Radius=min(R/2,C/2);
else
c0=(C+1)/2;
% Radius=min((R-1)/2,(C-1)/2);
end
Radius=min(r0-1,c0-1);

if NumberOfAngles==1
    Output=ACF(r0,c0:end);
else
ACF1=flipud(ACF(1:r0-1,c0:end));
ACF2=ACF(r0:end,c0:end);
ProfMat=zeros(NumberOfAngles*2,Radius);

for j=1:2
    if j==1
        y=ACF1';
    else
        y=ACF2;
    end
    
% CALCULATION OF ROTATIONAL MEAN
% Definition of angles
t=(pi/NumberOfAngles/2:pi/NumberOfAngles/2:pi/2);
   
% Matrix
y=y(1:Radius,1:Radius);
% Cycle between the 2nd and 2nd to last angles
[~, y1y]=size(y);

for i=1:NumberOfAngles
   rt=ceil(cos(t(i))*(1:Radius));
   ct=ceil(sin(t(i))*(1:Radius));
   profile=y((rt-1).*y1y+ct);

   if j==1
   ProfMat(NumberOfAngles+i,:)=profile;
   else
%    ProfMat(i,:)=profile;
   ProfMat(i,:)=[profile(2:end),profile(end)];  % excluding the central ACF point (0,0)
   end   
end

end


Output=[double(ACF(r0,c0)) sum(ProfMat)./(2*NumberOfAngles)];
% 
% OrientedProfiles=min_fw_Profile;
% OrientedProfiles(2,:)=max_fw_Profile;
% Angles=[min_th,max_th];

end


end

function bool=iseven(x)

if mod(x,2) == 0
bool=1;
else
bool=0;
end
end

function [param, fval, chisq]=simple_ICCS_Fit_ICS_1Dsingle(x,y,w0,Display, title1)
%fixed 

my=min(y);
My=max(y);
fun = @(Param) sum( (( (Param(1)+Param(2).*exp(-((x-0).*(x-0)./(Param(3)*Param(3)) ))) -y ).^2)./(abs(y))  ) ;
[param, chisqpar]=fminsearch( fun,[my My w0]);
param(3)=abs(param(3));
fval=(param(1)+param(2).*exp(-((x-0).*(x-0)./(param(3)*param(3)) )));
chisq=sum( (( (param(1)+param(2).*exp(-((x-0).*(x-0)./(param(3)*param(3)) ))) -y ).^2)./((param(2)^2)  )) ;

if Display==1
%     figure;
    plot(x,y,'o')
    hold on
    plot(x, fval, '--r')
    hold off
%     param(1)=param(1)./param(2);
    title(strcat(title1, '  w=',num2str(param(3),2),'   G0=',num2str(param(2),2) ));
end
end

function B=simpleICCS_Threshold(A,thr)
  
if length(thr)==1
B=A;
B(B<=thr)=0;
B(B>thr)=1;
else
B=A;
B(B>thr(2))=0;
B(B<=thr(1))=0;
B(B>0)=1;
end

end

function y=simpleICCS_smooth_simple(M,sm,n)
y=M;
if sm>0
filt = (1/(8+1/sm))*[1 1 1; 1 1/sm 1; 1 1 1]; % sm factor <=1 
    for i=1:n
    y = filter2(filt,y);
    end
end
    
end

