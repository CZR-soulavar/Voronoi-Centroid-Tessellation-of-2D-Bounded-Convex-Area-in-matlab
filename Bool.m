function [TX,stars]=Bool(farpre,faraft,X,A,B)
    vors=X{2,1};
    if X{3,1}==1 
        inf_posi=X{6,1};
        if inf_posi==1
            vors=[faraft;vors(2:end,:);farpre];
        elseif inf_posi==size(vors,1)
            vors=[vors(1:end-1,:);farpre;faraft];
        else
            vors=[vors(1:inf_posi-1,:);farpre;faraft;vors(inf_posi+1:end,:)];
        end
    end
        [posiNOx,posiNOy]=poly2cw(vors(:,1),vors(:,2));                    %因为polybool函数要求所使用的多边形顶点，要按照顺时针排序，poly2cw可以帮我们实现这个功能
        [X,Y]=polybool("intersection",A,B,posiNOx,posiNOy);                %将多边形与A,B组成的空间做布尔运算，得到交点，横坐标赋值给X，纵坐标赋值给Y
        polyin=polyshape(X,Y);                                             %polyin用于构造多边形变量                            
        [xpart,ypart]=centroid(polyin);                                    %计算布尔运算后多边形的质心
        stars=[xpart,ypart];
        TX=[transpose(X),transpose(Y)];
end