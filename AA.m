function[farpre,faraft]=AA(vorA,vorgen1,vorgen2)
    vora=vorA{4,1};
    gen=vorA{1,1};
    gen1=vorgen1{1,1};
    gen2=vorgen2{1,1};
    mid=zeros(4,2);
    mid(1,:)=(gen+gen1)/2;
    mid(2,:)=(gen+gen2)/2;
    mid(3,:)=2*vora-mid(1,:);
    mid(4,:)=2*vora-mid(2,:);
    relative1=mid(1,:)-vora;
    relative2=mid(2,:)-vora;
    relative3=mid(3,:)-vora;
    relative4=mid(4,:)-vora;
    angle=zeros(1,4);
    angle(1,1)=rad2deg(atan2(relative1(1,2),relative1(1,1)));
    angle(1,2)=rad2deg(atan2(relative2(1,2),relative2(1,1)));
    angle(1,3)=rad2deg(atan2(relative3(1,2),relative3(1,1)));
    angle(1,4)=rad2deg(atan2(relative4(1,2),relative4(1,1)));
    for i=1:1:4
        if angle(1,i)<0
            angle(1,i)=angle(1,i)+360;
        end
    end
    genrelative=gen-vora;
    anglegen=rad2deg(atan2(genrelative(1,2),genrelative(1,1)));
        if anglegen<0
            anglegen=anglegen+360;
        end
    [J,K]=sort(angle);
    t=4;
    for i=1:1:3
        if anglegen>J(1,i)&&anglegen<J(1,i+1)
            t=i;
        end
    end

    if t~=4
        pre=K(1,t);
        aft=K(1,t+1);
    elseif t==4
        pre=K(1,t);
        aft=K(1,1);
    end
    anglepre=deg2rad(angle(1,pre));
    angleaft=deg2rad(angle(1,aft));
    o=10^12;
    farpre=[o*cos(anglepre),o*sin(anglepre)]+vora;
    faraft=[o*cos(angleaft),o*sin(angleaft)]+vora;
end