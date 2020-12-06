%clear all,close all,clc
function [flag_i,flag_j]=core(I);
%I=imread('C:\Users\hp\Documents\MATLAB\corepoint\mine\db\12_2_2.bmp');
[m1 n1]=size(I)
me=mean2(I);
s1=std2(I);
dm=0;
ds=1;
flag_op=0;

%% normalisation 
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

%% Finding orientation field
if ~exist('orientsmoothsigma', 'var'), orientsmoothsigma = 0; end
    [rows,cols] = size(img);
        
    % Calculate image gradients.
    [Gx, Gy] = imgradientxy(img);
    
    %smoothing gradient
    if orientsmoothsigma
        sze = fix(6*orientsmoothsigma);   if ~mod(sze,2); sze = sze+1; end    
        f = fspecial('gaussian', sze, orientsmoothsigma);    
        Gx = filter2(f, Gx); % Smoothed sine and cosine of
        Gy = filter2(f, Gy); % doubled angles
    end
        
    % Coherence data for the image gradients     
    G2x = Gx.^2;       
    Gxy = Gx.*Gy;
    G2y = Gy.^2;
   
    %smoothing the coherence data by perform a weighted summation of the
    % data.
    blocksigma=5; %standard deviation of gaussian lowpass filter
    sze = fix(6*blocksigma);   if ~mod(sze,2); sze = sze+1; end    
    f = fspecial('gaussian', sze, blocksigma);
    G2x = filter2(f, G2x);
    Gxy = 2*filter2(f, Gxy);
    G2y = filter2(f, G2y);
   
    %sinusoidal component smoothing in orientation scale
    %computing sinusoidal components
    denom = sqrt(Gxy.^2 + (G2x - G2y).^2) + eps;
    sin2theta = Gxy./denom;            % Sine and cosine of doubled angles
    cos2theta = (G2x-G2y)./denom;
     
    %components convolution for orientation scale gaussian filter
    if orientsmoothsigma
        sze = fix(6*orientsmoothsigma);   if ~mod(sze,2); sze = sze+1; end    
        f = fspecial('gaussian', sze, orientsmoothsigma);    
        cos2theta = filter2(f, cos2theta); 
        sin2theta = filter2(f, sin2theta); 
    end
    %computing OF
    thee = pi/2 + atan2(sin2theta,cos2theta)/2;
    %figure;imshow(thee)
     
  %% finding edge leading to corepoint
    
    %orientation obtained generally forms v shaped images. The lowermost
    %point of the V shaped image coincides with the core of the fingerprint
    
    %Hence attempt is made to find the edge of the V shaped image and then
    %find lowermost point of the edge
    
       
    %converting the orientation field image to binary image
    %P is result of binarisation operation
    P=im2bw(thee);
    %figure;imshow(P)
        
    %finding edges of P using canny edge detection
    P=edge(P,'canny');
    %figure;imshow(P)
    
  %% Find lowest point of edge  
    %here P is the current name of the image that has result of canny edge
    %if any pixels valued 1 below (x,y) in 19*19 neighborhood, set P(x,y)=0
     
     %i and j are coordinates of current pixel, ie centre of
     %neighborhood window
     for i=10:m1-10
        for j=10:n1-10
            c=0;
        %value range of x is size of width of neighborhood window
          %since x ranges from i-9 to i+9 its size is: 
          %9(to left of i) + 1(current pixel i)+ 9(to right of i) =19
                  
        %value range of y is size of height of neighborhood window
          %since y ranges from j-9 to j+9 its size is: 
          %9(to left of j) + 1(current pixel i)+ 9(to right of j) =19
          
        %Hence size of window is 19*19
            for x=i-9:i+9        %value ranges from -9 to 9,ie 19
                for y=j-9:j+9    %value ranges from -9 tp 9, ie 19
                    if P(x,y)==1  
                        c=c+1;   %checks for number of pixels in neighborhood with value 1
                    end
                end
            end
            %if c=1, means that current pixel(centre of neighborhood
              %window) is the only 1 valued pixel in neighborhood
            %if c>1, means there exist other pixels with value 1 in
              %neighborhood. In this case we can set current pixel to 0
            if c>1 
                    P(i,j)=0;
            end
        end
     end
  %% Finding region of interest and selecting all relevant points within it
    
    %all the following operations doesn't effect image P
    %operations performed on J, which is binarised version of input image I
    J=I(:,:,1)>160;
    K=bwmorph(~J,'thin','inf');
    Kopen=imclose(K,strel('square',7));
    KopenClean=imfill(Kopen,'holes');
    KopenClean2=bwareaopen(KopenClean,7)
    KopenClean2([1 end],:)=0
    KopenClean2(:,[1 end])=0
    roi1=imerode(KopenClean2,strel('disk',10));
    %roi1 is the Region of interest obtained
    
  %% applying Region of Interest to select all relevant pixels
    %A is copy of roi1
    A=I;
    
    for i=1:m1
        for j= 1:n1
            if roi1(i,j)==0
                A(i,j)=0;
            else
                A(i,j)=1;
            end
        end
    end
    
    %if A value(ROI value)=0, the pixel lies out of ROI and need not be
    %considered
    for i=1:m1
        for j=1:n1
            if A(i,j)==0
                P(i,j)=0;  %lies outside ROI and hence set to 0
            end
        end
    end
    %figure;imshow(P)
 %% process to obtain corepoint from selected points in ROI
       
    %flag_i and flag_j are values of core coordinates finally passed to main function 
    flag_j=0;flag_i=0;
    
    maxd=300; 
    flag=0;
    co=0; 
    centerx=m1/2;centery=n1/2;
    
    %counting no of pixels with value 1 after applying ROI
    for i=1:m1
        for j=1:n1
            if P(i,j)==1
                co=co+1; %co is count
            end
        end
    end
    %if no of points obtained =1, only 1 pixel obtained, which is corepoint
 %% when co=1

    if co==1
        flag_op=1;
        for i=1:m1
            for j=1:n1
                if P(i,j)==1
                    flag_i=i;  %setting flag_i, flag_j to coordinates of pixel
                    flag_j=j;
                end
            end
        end
        
 %% when co=2
 
    %if no of points obtained =2, two pixels which maybe two corepoints or 
    %1 corepoint and 1 spurious value   
    elseif co==2
        flag_op=2;
        
        for i=1:rows
            for j=1:cols
                if P(i,j)==1 %only the two pixels pass if condition
                    
                    %centerx and centery are coordinates centroid of image
                    d=pdist2([i j],[centerx centery]); %distance between pixel and centroid
                    if d<maxd
                        maxd=d;
                        flagi=i;  %flagi and flagj set to coordinates of pixel that is closer to centroid
                        flagj=j;
                    end
                end
            end
        end
        %result of the previous for loop gives the coordinates of the pixel
        %among the two pixels under consideration which is closer to
        %centroid of image
        
        for i=1:rows
            for j=1:cols
                if (i~=flagi || j~=flagj) && P(i,j)==1 %if loop passes only the pixel that is further away from the centroid
                    d1=pdist2([i j],[flagi flagj]);    %to compute distance between the two pixels
                end 
            end
        end
                    if d1<30          %if both pixels within distance 30 of each other, taken as case of 2 corepoints of a whorl
                        for i=1:rows
                            for j=1:cols
                                if P(i,j)==1 
                                    f1_1=i; %loop repeats till f1_1 and f1_2 contains coords of bottom point
                                    f1_2=j;
                                    flag=1;
                                end
                            end
                        end
                        P(f1_1,f1_2)=0; %setting value of bottom point to 0
                                    
                    
                
                elseif d1>=30 %cases where distance>30, point closer to centroid taken as core 
                     
                    for i=1:rows
                        for j=1:cols
                            if i~=flagi && j~=flagj
                                P(i,j)=0; %setting point further away from centroid to 0
                            end
                        end
                    end
                end
            
        
 %% when co>=3
 
    %if no of points obtained >=3, when no of pixels selected after applying ROI>3.
    % the pixels are corepoint along with some spurious pixels
      
       elseif co>=3 %when co>=3
             
             f1_1=0;f1_2=0;
    %counting occurence of pixels in top and bottom half of image  
    flag_op=3;
 %figure;imshow(P)
 %to obtain lowermost 1 valued pixel in top half
 x1=co;
 f1_1=0;f1_2=0;flagi=0;flagj=0;
 maxd=300;
 
     while x1>1
         for i=2:rows-1
             for j=2:cols-1
                 if P(i,j)==1
                     d=pdist2([i j],[centerx centery]); %distance between pixel and centroid
                     if d<maxd
                        maxd=d;
                        flagi=i;  %flag_i and flag_j set to coordinates of pixel that is closer to centroid
                        flagj=j;
                     end
                 end
             end
         end
         flag=0;
         for i=2:rows-1
             for j=2:cols-1
                 if P(i,j)==1 && (i~=flagi || j~=flagj) && flag==0 %to select the second point. The first encountered point other than point nearest to centroid selected
                     f1_1=i;
                     f1_2=j; %f1_1 and f1_2 contain coordinates of second point
                     flag=1;
                 end
             end
         end
         d1=pdist2([f1_1 f1_2],[flagi flagj]); %finding distnace between both point       
         if d1<30
             if f1_1<flagi %if second point above point near centroid, setting point near centroid to 0
                 P(flagi,flagj)=0; 
             else    %when point near centoid above second point, set second point to 0
                 P(f1_1,f1_2)=0;
             end
         else %case where distance>30, set point further away to 0
             P(f1_1,f1_2)=0;
 
         end
         x1=x1-1; %decrement x1(param for while loop). Operation continues till only one point remains
     end
 
    end
 %% Store final values of flag_i and flag_j
 %after all operations only one pixel left with value 1
 %the coordinates of this pixel is selected as core
   for i=1:m1
        for j=1:n1
            if P(i,j)==1
                flag_i=i;
                flag_j=j;
            end
        end
  end
    
    if flag_i==0 && flag_j==0 %if no points detected, core set to centroid coordinate values
        flag_i=centerx;
        flag_j=centery;
    end
    %% printing final output
        
    %figure;imshow(I)
    %hold on
    %plot(flag_j,flag_i,'ro','linewidth',2);
    %op=insertShape(I,'circle',[flag_j flag_i 5 ],'Linewidth',3);
    
  %%
    
    