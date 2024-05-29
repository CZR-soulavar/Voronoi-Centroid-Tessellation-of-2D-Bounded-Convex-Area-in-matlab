function[mid,where]=BBB(A,B1,B2,p)
    gen=A{1,1};
    gen1=B1{1,1};
    gen2=B2{1,1};
    vorx=A{2,1};
    mid1=(gen+gen1)/2;
    mid2=(gen+gen2)/2;
    k1=(gen(1,2)-gen1(1,2))/(gen(1,1)-gen1(1,1));
    k2=(gen(1,2)-gen2(1,2))/(gen(1,1)-gen2(1,1));
    k1=-1/k1;
    k2=-1/k2;
    b1=mid1(1,2)-k1*mid1(1,1);
    b2=mid2(1,2)-k2*mid2(1,1);
    needdel=A{8,1};
    if p==1
        del=needdel(1,1);
    elseif p==2
        del=needdel(1,2);
    end
    vorx(del,:)=[999999,999999];
    x=vorx(:,1);
    y=vorx(:,2);
    yk1=k1*x+b1;
    yk2=k2*x+b2;
    deltay1=abs(yk1-y);
    deltay2=abs(yk2-y);
    t1=min(deltay1);
    t2=min(deltay2);
    if t1<t2
        genreal=gen2;
        where=B2(7,1);
    else 
        genreal=gen1;
        where=B1(7,1);
    end
    mid=(genreal+gen)/2;

end