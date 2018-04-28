%Inital Phase Start
%Finding the mean and standard deviation of all the 531 indices
load('20by531raw');
load('sfSeq.mat');
load('new');

load('mu4merby531.mat');
load('sigma4merby531.mat');
seq_mat=sfSeq(1:857,1:1);
%Aka=[0.05789;0.39521;0.05789];
%Initial Phase End
sms=sfSeq(1:580,1:1);
sz_sms=size(sms);
acsz_sms=sz_sms(1,1);
for z=1:acsz_sms
for Ai=1:531;
%Summarize Phase Start
%load('sequence_matrix_scop');
load('new');
keySet={'a','r','n','d','c','q','e','g','h','i','l','k','m','f','p','s','t','w','y','v'};
valueSet=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
mapObj = containers.Map(keySet,valueSet);

    proteinSeq=sms(z,:);
    chr = char(proteinSeq);
    L=length(chr);
    chr=chr';
    sum1=0;
    %count_Aka=1;
    new=new(:,Ai);
for i=1:L-2
    for j=i:i+2
        
        temp=chr(j);
        if(temp=='x' || temp=='X')
            break;
        end
        Ak= mapObj(temp);
        Lk= new(Ak);
        %wLk=Lk*Aka(count_Aka,1);
        %count_Aka=count_Aka+1;
        sum1=sum1+Lk;
         
    end
    avgsu=sum1/4;
    matrix(i,:)=avgsu;
    sum1=0;
    %count_Aka=1;
end
k(:,Ai)=matrix;

%Summarize Phase End

%Standardization Phase Starts
tempk=k(:,Ai);
numcols=size(tempk,2);
   for tempnorm=1:numcols
       seqOneByOne=tempk(:,tempnorm);
       sizeofOneCol=size(seqOneByOne);
       acSizeOfOneCol=sizeofOneCol(1:1);
       for tempAn=1:acSizeOfOneCol
		sso=seqOneByOne(tempAn,:);
        acNorm=(sso-mu4merby531(1,Ai))/sigma4merby531(1,Ai);
        matix_temp(tempAn,:)=acNorm;
       end
       An(:,tempnorm)=matix_temp;
   end
%Standardization Phase Ends

%Binning Phase Start
   %binsmatrix=hist(An,18);
    bins=[-4,-3.5,-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,1,2.5,3,3.5,4];
    bins=bins';
	Li=size(An);
    Li=Li(1:1);
	for idx = 1:numel(bins)
        if(bins(idx)==-4)
            count=sum(An<bins(idx));
        else
            count=sum((An<bins(idx)) & (An>=bins(idx-1)));
        end
         	prob=count/Li;
			binsmatrix(:,idx)=prob;
			if(idx==17)
			count=sum(An>bins(idx));
			prob=count/Li;
            binsmatrix(:,idx+1)=prob;
			end
        
    end
%Binning Phase End

if(Ai==1)
sfinalmatrix(Ai,:)=binsmatrix;
else
permanentma=horzcat(sfinalmatrix,binsmatrix);
sfinalmatrix=permanentma;
end
clear binsmatrix;
end
fm(z,:)=sfinalmatrix;
clear sfinalmatrix;
clear An;
clear k;
end