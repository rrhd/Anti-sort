clc, clear all
S=[-26, -97, 49, -61, -43, 6, -90, -71, 33, -69, -61, 15, 51, -21, -5, -90, -78, 68, 79, 27, -43, -34, 20, -36, 96, 29, -65, 88, -17, 58, -70, -47, -34, 86, -35, 39, 76, 0, -11, -66, -86, 53, -1, 17, 51, 35, 90, -23, 50, -28, 49, 59, -92, -48, -98, -3, -32, 24, 30, -28, 64, 87, 9, -7, 79, -85, 56, -50, 30, 72, 94, 67, -30, 4, -27, -52, 87, 96, -68, 57, 28, 36, 4, -37, -58, 84, 79, 89, 73, 26, 8, -92, 80, 64, -13, 77, 4, -64, -6, 35];
s=[6,98,88,43,34,76,89,11,62,10,58,33,72,68,91,97,81,55,18,17,47,75,42,22,94,29,51,60,63,39,30,2,70,83,37,66,25,5,45,41,82,38,65,50,100,32,23,16,67,8,49,35,26,56,87,79,71,15,12,48,9,84,13,54,78,92,77,14,44,53,46,99,52,85,93,4,69,7,86,24,59,31,19,64,96,1,61,21,28,95,3,57,90,73,20,27,36,40,80,74];
F = @(S,s) sum(abs(diff(S(s))));
cost = F(S,s);
T = 1;
T_min = 0.000001;
alpha=0.9999;
while T > T_min
    i=1;
    while i<=100
        k=randperm(100);
        sNew=s;
        sNew(k(1))=s(k(2));
        sNew(k(2))=s(k(1));
        newCost = F(S,sNew);
        ap = exp((newCost-cost)/T);
        if ap > rand
            s = sNew;
            cost = newCost;
        end
        i = i+1;
    end
    T = T*alpha;
end

F = @(S,s) sum(abs(diff(S(s))));
S=[-26, -97, 49, -61, -43, 6, -90, -71, 33, -69, -61, 15, 51, -21, -5, -90, -78, 68, 79, 27, -43, -34, 20, -36, 96, 29, -65, 88, -17, 58, -70, -47, -34, 86, -35, 39, 76, 0, -11, -66, -86, 53, -1, 17, 51, 35, 90, -23, 50, -28, 49, 59, -92, -48, -98, -3, -32, 24, 30, -28, 64, 87, 9, -7, 79, -85, 56, -50, 30, 72, 94, 67, -30, 4, -27, -52, 87, 96, -68, 57, 28, 36, 4, -37, -58, 84, 79, 89, 73, 26, 8, -92, 80, 64, -13, 77, 4, -64, -6, 35];
cost=F(S,1:100);
costNew=0;
K=S;
tol=1;
iter=1;
while abs(costNew-cost)>tol
    if iter>1
        K=newK;
        cost=costNew;
    end
    B=K(iter:end);
    DiffMat=zeros(length(B),length(B));
    for i=1:length(B)
        for j=1:length(B)
            DiffMat(i,j)=abs(B(j)-B(i));
        end
    end
    DiffMat=DiffMat - diag(diag(DiffMat));
    DiffMat=tril(DiffMat);
    [C,Ind1]=max(DiffMat);
    maxDiff=max(C);
    Ind2=find(C==maxDiff);
    MaxDiffMax=zeros(length(Ind2),1);
    for i=1:length(Ind2)
        MaxDiffMax(i)=sum(DiffMat(Ind1(Ind2(i)),:))+sum(DiffMat(:,Ind2(i)));
    end
    [C,Ind3]=max(MaxDiffMax);
    if length(find(C==MaxDiffMax))>1
        keyboard
    end
    sum1=sum(DiffMat(Ind1(Ind2(Ind3)),:));
    sum2=sum(DiffMat(:,Ind2(Ind3)));
    newK=K;
    if sum1>=sum2
        if Ind2(Ind3)<Ind1(Ind2(Ind3))
            newK(Ind2(Ind3))=[];
            newK=[newK(1:iter-1),K(Ind2(Ind3)),newK(iter:end)];
            costNew=F(newK,1:length(newK));
        else
        end
    else
        
        if Ind2(Ind3)<Ind1(Ind2(Ind3))
            newK(Ind1(Ind2(Ind3)))=[];
            newK=[newK(1:iter-1),K(Ind1(Ind2(Ind3))),newK(iter+1:end)];
            costNew=F(newK,1:length(newK));
        else
        end
    end
    iter=iter+1;
end

