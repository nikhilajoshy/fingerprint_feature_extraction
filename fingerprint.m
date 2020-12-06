function C=fingerprint(I)
%{
clear all,close all,clc
%p=input();
I=imread('C:\Users\hp\Documents\MATLAB\corepoint\mine\db\1_2_2.bmp');
%figure;imshow(I)
%}
%% cropping image based on corepoint as centre
[Cx,Cy]=core(I); %Cx and Cy are coordinates of core 
I=histeq(I);

xmin=Cx-50;
ymin=Cy-50;
I=imcrop(I,[ymin,xmin,100,100]); %cropping I to image of dimension 100*100 
%figure;imshow(I)


%% normalisation
me=mean2(I);
s1=std2(I);
dm=0;
ds=1;
[m1,n1]=size(I); 
 for i=1:m1
     for j=1:n1
         if I(i,j)>me
             
         t= dm+sqrt(double((ds*((I(i,j)-me)^2))/s1));
         else
             t= dm-sqrt(double((ds*((I(i,j)-me)^2))/s1));
         end
         img(i,j)=t;
     end
 end
%%  

%image enhancement
I=histeq(img);
I=imadjust(I);

%binarisation
J=im2bw(I,0.5);

%thinning
K=bwmorph(~J,'thin','inf');
%figure;imshow(K)


%applying minutiae on image to find no: of neighboring pixels
%minutiae is function that operates on 3*3 windows, returning no:of
%neighbours of center pixel
fun=@minutie;
L = nlfilter(K,[3 3],fun); %L is no: of pixels in neighborhood of center pixel

%finding centroid of termination
LTerm=(L==1); %if L=1, pixel is a ridge termination
LTermLab=bwlabel(LTerm);
propTerm=regionprops(LTermLab,'Centroid');
CentroidTerm=round(cat(1,propTerm(:).Centroid)); %concatenating pixel to array of Termination points

%centroid of bifurcation
LBif=(L==3); %if L=3, pixel is a ridge bifurcation
LBifLab=bwlabel(LBif);
propBif=regionprops(LBifLab,'Centroid','Image');
CentroidBif=round(cat(1,propBif(:).Centroid)); %concatenating pixel to array of Bifurcation points

D=4;

Distance=DistEuclidian(CentroidBif,CentroidTerm); %finding distances between pairs of terminations and bifurcations point
SpuriousMinutae=Distance<D;  %if distance<D
[i,j]=find(SpuriousMinutae); %find index values of termination and bifrcation point in respective arrays
CentroidBif(i,:)=[];         %remove termination and bifurcation points from respective arrays  
CentroidTerm(j,:)=[];

Distance=DistEuclidian(CentroidBif); %finding distances between all pairs of bifurcations
SpuriousMinutae=Distance<D;
[i,j]=find(SpuriousMinutae);
CentroidBif(i,:)=[];

Distance=DistEuclidian(CentroidTerm); %finding distances between all pairs of termination
SpuriousMinutae=Distance<D;
[i,j]=find(SpuriousMinutae);
CentroidTerm(i,:)=[];



%%finding Region of interest
v=0;
Kopen=imclose(K,strel('square',7)); %operations to find ROI includes:
KopenClean= imfill(Kopen,'holes');   %closing,filling, area open, erosion
KopenClean=bwareaopen(KopenClean,5);
KopenClean([1 end],:)=0;
KopenClean(:,[1 end])=0;
ROI=imerode(KopenClean,strel('disk',10)); %resulting region of interest
%figure;imshow(ROI)
%% suppressing minutiae to ROI so only minutiae with ROI considered

%in this method a new array Z is created where only points corresponding to
%termination points have value 1. It is then multiplied with ROI which has
%value 1 in only region of interest. Result stored in ZTerm thus only
%contains termination points in region of interest. Similar operation is
%performed for bifurcation points

[m,n]=size(I(:,:,1));
indTerm=sub2ind([n,m],CentroidTerm(:,1),CentroidTerm(:,2)); %indTerm is linear index equivalent to  row and col subscripts of CentroidTerm 
Z=zeros(n,m); %creating new array of size n,m 
Z(indTerm)=1; %setting values corresponding to index position of minutiae points to 1
ZTerm=Z.*ROI'; %multiplying by ROI, so that only termination points(value 1) in ROI(ROI has value 1 in region of interest) will have final val 1
[CentroidTermX,CentroidTermY]=find(ZTerm); %x, y coordinates of array ZTerm gives coordinates of termination points in ROI

indBif=sub2ind([n,m],CentroidBif(:,1),CentroidBif(:,2)); %indBif is linear index of bifurcation points
Z=zeros(n,m);
Z(indBif)=1;
ZBif=Z.*ROI';
[CentroidBifX,CentroidBifY]=find(ZBif); %x and y coordinates of ZBif gives coordinates of Bifurcation points in ROI

%TO PRINT TERMINATION AND BIFURCATION POINTS
%{
figure;imshow(K)
set(gcf,'position',[1 1 1000 1000]);
hold on
plot(CentroidTermX,CentroidTermY,'ro','linewidth',2)
plot(CentroidBifX,CentroidBifY,'go','linewidth',2)
%}
%% Storing x and y coordinates of minutiae points to arrays MX and MY
%MX contains all x coordinates of terminations followed by x coordinates of
%bifurcations.
%Similar method for y coordinates where MY is final set of Y coordinates

s1=size(CentroidTermX);
s2=size(CentroidBifX);

%s1 no of terminations
%s2 no of bifurcationss
s1=s1(1);
s2=s2(1);

%MX is set of X coordinates of all minutiae points detected
%MX(1) and MY(1) correspond to coordinates of first minutiae point
if s1>0
    for i=1:s1
        MX(i)=CentroidTermX(i); %first storing all termination x coordinates
        MY(i)=CentroidTermY(i); %first storing all termination y coordinates
    end
    i=i+1;
else
    i=1;
end
v=1;
x=i;
if s2>0
    for j=i:i+s2-1
        MX(j)=CentroidBifX(v); %concatenating values of bifurcation x coordinates to existing termination x coordinates 
        MY(j)=CentroidBifY(v); %concatenating values of bifurcation y coordinates to existing termination y coordinates
        v=v+1;
    end
end

%%
%Storing all coordinate values as binary string C
C='';
for i=1:j %the for loop concatenates binary string of x and y coordinates of ith minutiae at each iteration
    %A contains 7 bit binary version of x coordinate
    A=dec2bin(MX(i),7);
    for x=1:7
        A(x)=num2str(A(x)); %converting binary to string
        C=strcat(C,A(x)); %x coordinate of ith coordinate concatenated to C
    end
    %B contains 7 bit binary version of y coordinate
    B=dec2bin(MY(i),7);
    for x=1:7
        B(x)=num2str(B(x)); %converting binary to string
        C=strcat(C,B(x)); %y coordinate of ith coordinate concatenated to C
    end
end %C now contains the coordinates of all minutiae points and is passed to main function for watermarking