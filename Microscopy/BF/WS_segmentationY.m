% Watershedding by Denis Tsygankov

function G=WS_segmentationY(F,Hcr)
    tic
    [SF, SI]=sort(F(:));
    z=find(SF>0,1,'first');
    MY=size(F,1);
    MX=size(F,2);
    L=MX*MY;
    G=zeros(size(F));
    PV=zeros(L,1);

    cnt=1; 
    
    x=ceil(SI(L)/MY);
    y=SI(L)-MY*(x-1);
    G(y,x)=cnt;
    PV(cnt)=F(y,x);
    
    for i=(L-1):(-1):z
        x=ceil(SI(i)/MY);
        y=SI(i)-MY*(x-1);
    
        A=G((y-1):(y+1),(x-1):(x+1));
        B=unique(A(:)); LB=length(B);
     
        if LB==1
            cnt=cnt+1;
            G(y,x)=cnt;
            PV(cnt)=F(y,x);                
        elseif LB==2
            G(y,x)=B(2);                
        elseif LB>=3
            bV=F(y,x);
            tV=PV(B(2:LB));
            [Hv,Hi]=sort(tV-bV);
            if Hv(LB-2)<Hcr
                for k=1:(LB-2)
                    G(G==B(Hi(k)+1))=B(Hi(LB-1)+1);                    
                end
                G(y,x)=B(Hi(LB-1)+1);
            end
        end
    end
%     toc




