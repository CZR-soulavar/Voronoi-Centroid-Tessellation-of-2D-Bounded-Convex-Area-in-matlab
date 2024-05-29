%% 函数输入区，输出区，输入的值有：边界的横纵坐标AB、初始值点XY，输出的值有Bright
function[newgenspre,CHECK]=Vofs(A,B,X,Y,dx)
%tic
    N=size(X,1);
    Z=[X,Y];
    newgens=zeros(N,2)+9999;                                                   % 用newgens记录所有新的生成点的坐标
    newgenspre=zeros(N,2);

 while sum(abs(newgens(:,1)-newgenspre(:,1)))>dx
        zhishi=sum(abs(newgens(:,1)-newgenspre(:,1)));
%% 变量与初始值定义区
        [C,V]=voronoin(Z);                                                         % 使用生成点进行维诺划分，其中C是维诺细胞顶点的坐标，V是维诺细胞顶点的坐标索引
        Infs=transpose(find(cellfun(@(x)any(x == 1),V)));                          % 通过维诺细胞顶点的坐标索引V，找到所有无界维诺细胞在V中所处的位置，并代入Infs
        xulie=1:1:N;                                                               % 获取序列1到N，用其删除无限维诺细胞的位置，剩下的全是有限维诺细胞
        Fins=setdiff(xulie,Infs);                                                  % 获取有限维诺细胞在V中的位置，并存入Fins
        ENDS=cell(8,N);                                                            
%为函数设置元胞数组变量，用于存储每个细胞的信息
%第一行为生成点坐标
%第二行为顶点坐标
%第三行为1代表该细胞为是无界维诺细胞，为0代表是有界维诺细胞。
%当第三行为0时，第四行为0；当第三行为1时，第四行为无限顶点前后顶点的坐标
%当第三行为0时，第五行为0；当第三行为1时，若该细胞仅有一个有限顶点，则第五行为1，若该细胞有多个有限顶点，则第五行为2
%当第三行为0时，第六行为0；当第三行为1时，会返回该细胞无限顶点所在的位置

        vorpre=zeros(N,2);                                                         % 用vorpre记录Inf_count个无界维诺细胞的无限顶点旁的前顶点
        voraft=zeros(N,2);                                                         % 用voraft记录Inf_count个无界维诺细胞的无限顶点旁的后顶点
        vornumber=zeros(1,N);                                                      % 用vornumber记录Inf_count个无界维诺细胞的顶点个数
        atype=[];                                                                  % 用atype记录所有只有一个有限顶点的无界维诺细胞的编号
        newgenspre=newgens;
        newgens=zeros(N,2);                                                        % 用newgens记录所有新的生成点的坐标
        CHECK=cell(1,N);                                                           % 用CHECK元胞数组检查
%% ENDS补完区，按照上述定义变量补完ENDS这一元胞数组，方便后续进行类型判断，数据取用
    for i=1:1:N                                                                % 对于第i个维诺细胞来说
            ENDS(1,i)={Z(i,:)};                                                    % ENDS的第一行的第i个元素，是第i个维诺细胞生成点的坐标，来源于Z
            ENDS(7,i)={i};                                                         % ENDS的第七行的第i个元素是第i个维诺细胞的编号
            vor=[];                                                                % 将vor空矩阵以储存第i个维诺细胞的顶点坐标
            vorindex=cell2mat(V(i,1));                                             % 取出第i个细胞在C中对应的顶点位置

            for j=1:1:length(vorindex)                                             % 对于第i个细胞的第j个顶点来说
                vor=[vor;C(vorindex(1,j),:)];                                      % 将第i个细胞的第j个顶点的坐标放在vor数组的第j行
            end                             

            ENDS(2,i)={vor};                                                       % ENDS的第二行的第i个元素是第i个维诺细胞的所有顶点坐标

            if ismember(i,Infs)                                                    % 如果第i个细胞是无界维诺细胞
                ENDS(3,i)={1};                                                     % 则将ENDS的第三行第i个元素置为1
            else                                                                   % 如果第i个细胞是有界维诺细胞
                ENDS(3,i)={0};                                                     % 则将ENDS的第三行第i个元素置为0                                             
            end
            vornumber(1,i)=size(vor,1);                                            % 得到第i个维诺细胞的顶点个数为vornumber
            ENDS(4,i)={0};                                                         % 为ENDS第四行提前预设值，防止空元胞数组出现
            ENDS(5,i)={0};                                                         % 为ENDS第五行提前预设值，防止空元胞数组出现
            ENDS(6,i)={0};                                                         % 为ENDS第六行提前预设值，防止空元胞数组出现
            ENDS(8,i)={0};                                                         % 为ENDS第八行提前预设值，防止空元胞数组出现

            if ismember(i,Infs)                                                    % 如果第i个细胞是无界维诺细胞
                vornumber(1,i)=size(vor,1);                                        % 得到第i个无界维诺细胞的顶点个数为vornumber
                Inf_posi=find(vorindex==1);                                        % 找到第i个无界维诺细胞的无穷顶点在索引中的位置，代入Inf_posi
                ENDS(6,i)={Inf_posi};                                              % 将ENDS第六行的第i个元素为维诺细胞无限顶点所在位置
                if Inf_posi==1                                                     % 如果无限顶点在第一个位置，则前顶点和后顶点分别是该无界维诺细胞顶点序列的倒数第一个和第二个
                    vorpre(i,:)=C(vorindex(1,end),:);                                            %vorpre=该细胞的倒数第一个顶点坐标
                    voraft(i,:)=C(vorindex(1,2),:);                                              %voraft=该细胞的第一个顶点坐标
                    ENDS(8,i)={[vornumber(1,i),2]};                                              %ENDS第八行的第i个元素为维诺细胞无限顶点前后顶点所在位置分别为：顶点数组长度，2
                elseif Inf_posi==length(vor)                                       % 如果无限顶点在最后的位置，则前顶点和后顶点分别是该无界维诺细胞序列的倒数第二个和第一个
                    vorpre(i,:)=C(vorindex(1,end-1),:);                                          %vorpre=该细胞的倒数第二个顶点坐标
                    voraft(i,:)=C(vorindex(1,1),:);                                              %voraft=该细胞的第一个顶点坐标
                    ENDS(8,i)={[vornumber(1,i)-1,1]};                                            %ENDS第八行的第i个元素为维诺细胞无限顶点前后顶点所在位置分别为：顶点数组长度-1，1
                else                                                               % 如果无限顶点不在第一个也不在最后一个，则前顶点和后顶点分别是该无限顶点的前一个和后一个
                    vorpre(i,:)=C(vorindex(1,Inf_posi-1),:);                                     %vorpre=该细胞无限顶点的前顶点的坐标
                    voraft(i,:)=C(vorindex(1,Inf_posi+1),:);                                     %voraft=该细胞无限顶点的后顶点的坐标
                    ENDS(8,i)={[Inf_posi-1,Inf_posi+1]};
                end

                if vornumber(1,i)==2                                               % 如果第i个维诺细胞是无界维诺细胞且仅有2个顶点（此时已经是在无界维诺细胞的判断条件下了，所以用两个判断没问题）
                    atype=[atype,i];                                               % 将此维诺细胞的编号代入atype
                    ENDS(5,i)={1};                                                 % 将ENDS第五行的第i个元素置1，说明这是A类无界维诺细胞
                    ENDS(4,i)={[vorpre(i,:)]};                                     % 将ENDS第四行的第i个元素置为无限顶点前顶点的坐标（因为只有一个有限顶点，前后坐标一样）
                else                                                               % 如果第i个维诺细胞是无界维诺细胞且有>2个顶点
                    ENDS(5,i)={2};
                    ENDS(4,i)={[vorpre(i,:);voraft(i,:)]};                         % 将ENDS第四行的第i个元素置为无限顶点前后顶点的坐标（分为两列）
                end
    end
   
 end 

%% 边界细胞处理区
        for i=1:1:N
          if ENDS{5,i}~=0                                                        % 当该细胞属于A或B型维诺细胞时
                distancespre=pdist2(vorpre(i,:),Z);                                % 计算所有生成点和无限顶点前顶点的距离d1
                distancesaft=pdist2(voraft(i,:),Z);                                % 计算所有生成点和无限顶点后顶点的距离d2
                distvorpre=pdist2(vorpre(i,:),ENDS{1,i});                          % 计算无限顶点前顶点和该细胞生成点的距离d3
                distvoraft=pdist2(voraft(i,:),ENDS{1,i});                          % 计算无限顶点后顶点和该细胞生成点的距离d4
                distrelativepre=abs(distancespre-distvorpre);                      % 计算|d1-d3|的值并存入distrelativepre中
                distrelativeaft=abs(distancesaft-distvoraft);                      % 计算|d2-d4|的值并存入distrelativeaft中
                [least3pre,least3_posipre]=sort(distrelativepre);                  % 对计算出的|d1-d3|值进行排序，最小的3个放在least3_posipre中，其中least3_pre是最小三个在|d1-d3|中的索引
                [least3aft,least3_posiaft]=sort(distrelativeaft);                  % 对计算出的|d2-d4|值进行排序，最小的3个放在least3_posiaft中，其中least3_aft是最小三个在|d2-d4|中的索引
                least3_posipre=least3_posipre(1:3);                                % |将无限顶点前顶点到其该细胞生成点距离-无限顶点前顶点到其他细胞生成点距离|的三个点索引取出放在least3_posipre中
                least3_posiaft=least3_posiaft(1:3);                                % |将无限顶点后顶点到其该细胞生成点距离-无限顶点后顶点到其他细胞生成点距离|的三个点索引取出放在least3_posiaft中
                least2_posipre=setdiff(least3_posipre,i);                          % 将该细胞生成点从三个点中排除，剩下的两个就是该顶点的其他两个生成点的索引
                least2_posiaft=setdiff(least3_posiaft,i);                          % 将该细胞生成点从三个点中排除，剩下的两个就是该顶点的其他两个生成点的索引
                vorgenpre1=least2_posipre(1,1);                                    % 将该细胞的前顶点的两个非本细胞生成点之一的索引代入vorgenpre1
                vorgenpre2=least2_posipre(1,2);                                    % 将该细胞的前顶点的两个非本细胞生成点之一的索引代入vorgenpre2
                vorgenaft1=least2_posiaft(1,1);                                    % 将该细胞的后顶点的两个非本细胞生成点之一的索引代入vorgenaft1
                vorgenaft2=least2_posiaft(1,2);                                    % 将该细胞的后顶点的两个非本细胞生成点之一的索引代入vorgenaft2

            if  ENDS{5,i}==2                                               % 如果第i个细胞是B型维诺细胞

                if isempty(intersect(least2_posipre,atype))==0                     % 对于无限顶点前顶点来说，如果其他两个生成点所属细胞中有一个属于A型维诺细胞，则必然该顶点是和该维诺细胞生成的--BAB型顶点
                     wherepre=intersect(least2_posipre,atype);                     % 对于BAB型中的B类细胞的无限顶点前顶点来说，找出那个A型维诺细胞的位置
                     [midpre]=BAB(ENDS(:,i),ENDS(:,wherepre));                     % 使用函数BAB代入第i个维诺细胞的全部信息与找到的A型维诺细胞的全部信息，得到前顶点对应的中点
                elseif isempty(intersect(least2_posipre,Fins))==0                  % 对于无限顶点前顶点来说，如果其他两个生成点所属细胞中有一个属于C型维诺细胞，则必然该顶点是和另外一个生成点所属维诺细胞生成的---BBC型顶点
                     wherepre=setdiff(least2_posipre,intersect(least2_posipre,Fins));%对于BBC型中的B类细胞的无限顶点前顶点来说，找出那个C型维诺细胞的位置，并取另一个B作为生成点
                     [midpre]=BBC(ENDS(:,i),ENDS(:,wherepre));                     % 使用函数BBC代入第i个维诺细胞的全部信息与找到的C型维诺细胞的全部信息，得到前顶点对应的中点
                else                                                               % 对于无限顶点前顶点来说，如果不属于上述的任何一种情况，则需要再议，但这种情况在大数量级下存在的可能性微乎其微---BBB型顶点
                    p=1;                                                           % p=1，代表此时使用的顶点是第i个维诺细胞的无限顶点前顶点
                    [midpre,wherepre]=BBB(ENDS(:,i),ENDS(:,vorgenpre1),ENDS(:,vorgenpre2),p);% 对于BBB型中的B类细胞的无限顶点前顶点来说，需要用该细胞生成点与另外两个生成点做垂直平分线，并确定是否有顶点处于该垂直平分线上，如果在，则该顶点生成点是另外一个
                end     

                if isempty(intersect(least2_posiaft,atype))==0                     % 对于无限顶点后顶点来说，如果其他两个生成点所属细胞中有一个属于A型维诺细胞，则必然该顶点是和该维诺细胞生成的--BAB型顶点
                     whereaft=intersect(least2_posiaft,atype);                     % 对于BAB型中的B类细胞的无限顶点后顶点来说，找出那个A型维诺细胞的位置
                     [midaft]=BAB(ENDS(:,i),ENDS(:,whereaft));                     % 使用函数BAB代入第i个维诺细胞的全部信息与找到的A型维诺细胞的全部信息，得到后顶点对应的中点
                elseif isempty(intersect(least2_posiaft,Fins))==0                  % 对于无限顶点后顶点来说，如果其他两个生成点所属细胞中有一个属于C型维诺细胞，则必然该顶点是和另外一个生成点所属维诺细胞生成的---BBC型顶点
                     whereaft=setdiff(least2_posiaft,intersect(least2_posiaft,Fins)); %对于BBC型中的B类细胞的无限顶点后顶点来说，找出那个C型维诺细胞的位置，并取另一个B作为生成点
                     [midaft]=BBC(ENDS(:,i),ENDS(:,whereaft));                     % 使用函数BBC代入第i个维诺细胞的全部信息与找到的C型维诺细胞的全部信息，得到后顶点对应的中点
                else                                                               % 对于无限顶点后顶点来说，如果不属于上述的任何一种情况，则需要再议，但这种情况在大数量级下存在的可能性微乎其微---BBB型顶点
                    p=2;                                                           % p=2，代表此时使用的顶点式第i个维诺细胞的无限顶点后顶点
                    [midaft,whereaft]=BBB(ENDS(:,i),ENDS(:,vorgenaft1),ENDS(:,vorgenaft2),p);% 对于BBB型中的B类细胞的无限顶点后顶点来说，需要用该细胞生成点与另外两个生成点做垂直平分线，并确定是否有顶点处于该垂直平分线上，如果在，则该顶点生成点是另外一个
                end     

                     [farpre,faraft]=Calculation(ENDS(:,i),midpre,midaft);         % 将前后中点代入，分别与前后顶点组成直线相交后得到交点，后反向延长得到远处点farpre和faraft

           elseif ENDS{5,i}==1                                                     % 当该细胞属于A型维诺细胞时，不论该细胞的前后顶点生成点属于AAA还是BAB型，都是要使用极坐标方法处理
                     [farpre,faraft]=AA(ENDS(:,i),ENDS(:,vorgenpre1),ENDS(:,vorgenpre2));%使用极坐标方法处理A型维诺细胞，得到点
            end

         else
                     farpre=[0,0];                                                 % 如果该维诺细胞是有限维诺细胞，则远处前顶点设为[0,0]
                     faraft=[0,0];                                                 % 如果该维诺细胞是有限维诺细胞，则远处后顶点设为[0,0]
          end
            [TX,stars]=Bool(farpre,faraft,ENDS(:,i),A,B);                          % 用Bool函数进行布尔运算，输出的star是第i个维诺细胞交区域的质心（不论有界无界都能进行计算），TX是第i个维诺细胞的新顶点
            newgens(i,:)=stars;                                                    % 用newgens的第i行记录N个维诺细胞新交区域的质心
            CHECK(1,i)={TX};                                                       % 用CHECK元胞数组的第i行记录N个维诺细胞交区域的顶点          
        end
        X=newgens(:,1);                                                            % 将新生成点的横坐标赋值给X
        Y=newgens(:,2);                                                            % 将新生成点的纵坐标赋值给Y
        Z=[X,Y];                                                                   % 将[X,Y]代入Z
 end
        %voronoi(X,Y);                                                              % 为最后维诺质心镶嵌的结果绘制维诺图
        %hold on
        %plot(A,B);
        %axis equal                                                                 % 横纵等距离    
        %toc
end