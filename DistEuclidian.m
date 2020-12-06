function D=DistEuclidian(dataset1,dataset2) %dataset is the minutiae array passed
switch nargin
    case 1 %where datasets are either 2 copies of same terminations array or bifurcation array
        [m1,n1]=size(dataset1);
        m2=m1; %since both data sets are the same
        D=zeros(m1,m2);
        for i=1:m1
            
            for j=1:m2
                if i==j
                    D(i,j)=NaN;
                else %computing distance between ith coordinate in dataset to jth coordinate
                    D(i,j)=sqrt((dataset1(i,1)-dataset1(j,1))^2+(dataset1(i,2)-dataset1(j,2))^2);
                end
            end
        end
    case 2 %case where 1 dataset is termination array and the other is bifurcation array
        [m1,n1]=size(dataset1);
        [m2,~]=size(dataset2);
        D=zeros(m1,m2);
        for i=1:m1
            for j=1:m2 %computing distance between ith coordinate in termination array and jth coordinate in bifurcation array
                D(i,j)=sqrt((dataset1(i,1)-dataset2(j,1))^2+(dataset1(i,2)-dataset2(j,2))^2);
            end
        end
    otherwise %default
        error('only one or two input arguments')
end

