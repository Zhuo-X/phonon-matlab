NA=6.02214076*10^23;
a = 1;   
MB=10.8*1/(NA)/1000;%*11/10
MN=14.1/(NA)/1000;%*11/10  #mass kf N(kg)
MB1=11*1/(NA)/1000;
MN1=15*1/(NA)/1000
a1=(1/MB1-1/MB)*MB;
a2=(1/sqrt(MN)*1/sqrt(MB1)-1/sqrt(MN)*1/sqrt(MB))*sqrt(MN)*sqrt(MB);
c=0.5;%11的浓度
g1=c*(1-11/(11*c+10*(1-c)))+(1-c)*(1-10/(11*c+10*(1-c)));
g2=c*(1-11/(11*c+10*(1-c)))^2+(1-c)*(1-10/(11*c+10*(1-c)))^2;
fn=diag([352.446813, 80.26905, 84.19438,54.97067986,-29.0154186, -6.69709,-41.07980856, 37.63512882,7.033555,7.674426014,1.233675998,7.09764]);
fb=fn;
[FA,FB] = rotate(2/3*pi, [1;0;0]*a);
[SA, SB] = rotate(1/3*pi, [3/2; sqrt(3)/2; 0]*a);
[TA, TB] = rotate(2/3*pi, [1; sqrt(3); 0]*a) ;        % 第三近邻
[LA1, LB1] = rotate(2/3*pi, [2.5; sqrt(3)/2; 0]*a);   % 第四近邻，两种角度
[LA2, LB2] = rotate(2/3*pi, [2.5; -sqrt(3)/2; 0]*a);
LA =cat(3,LA1,LA2) ;                                             % 第四近邻两种情况进行合并，可省略
LB = cat(3,LB1,LB2);
[KAB1,KBA1] = rotate(2/3*pi, fb(1:3, 1:3));%每个原子的K
[KAA,KBB] = rotate(1/3*pi, K(1/6*pi, fb(4:6, 4:6)));
[KAB3,KBA3] = rotate(2/3*pi, K(1/3*pi, fb(7:9, 7:9)));
[KAB4f,KBA4f] = rotate(2/3*pi, K(acos(2.5/sqrt(7)), fb(10:12, 10:12)));
[KAB4s,KBA4s] = rotate(2/3*pi, K(2*pi-acos(2.5/sqrt(7)), fb(10:12, 10:12)));
KAB4 = cat(3,KAB4f,KAB4s);
KBA4 = cat(3,KBA4f,KBA4s);
n = 10 ;        % 沿第一布里渊区的最小重复单元的边界每条边单位长度取100个k点进行计算（主要是为了使点均匀分布，同样可以每条边取n 个点进行计算）
result1 = zeros([(30 + 17) * n, 7]) ;  % 解的矩阵
result = zeros([(30 + 17) * n, 7]);
result2 = zeros([(30 + 17) * n, 7]);
result3 = zeros([(30 + 17) * n, 7]);
result4= zeros([(30 + 17) * n, 7]);
result5= zeros([(30 + 17) * n, 7]);
result6= zeros([(30 + 17) * n, 7]);
for i=1:1:((30 + 17) * n) 
    if i < n * 17                       
        kx = i*2*pi/3/a/(n * 17);
        ky = 0;
    elseif i < (10 + 17) * n
        kx = 2 * pi / 3 / a;
        ky = (i-n * 17) / (10 * n-1) * 2 * pi / 3 / a/sqrt(3);
    else
        kx = 2 * pi / 3 / a - (i-(10 + 17) * n)/(n * 20-1) * 2 * pi / 3 / a;
        ky = kx/sqrt(3);
    end
    k = [kx, ky, 0];
    D = zeros(6);
    D1 = zeros(6);
    D2=zeros(6);
    D3=zeros(6);
    D4=zeros(6);
    D5=zeros(6);
    DAAs = 0;
    DBBs = 0;
    DBAs = 0;
    DABs = 0;
    for j = 1:3
        DAAs = DAAs+KAA(:,:,j)*exp(1j*(k* -SA(:,:,j))) + KAA(:,:,j+3)*exp(1j*(k* -SA(:,:,j+3)));
        DBBs = DBBs + KBB(:,:,j) * exp(1j * (k* -SB(:,:,j))) + KBB(:,:,j+3)*exp(1j*(k* -SB(:,:,j+3)));
        DABs = (DABs + KAB1(:,:,j) * exp(1j * (k* -FA(:,:,j))) + KAB3(:,:,j) * exp(1j * (k* -TA(:,:,j))) + ...
                KAB4(:,:,j) * exp(1j * (k* -LA(:,:,j))) + KAB4(:,:,j+3) * exp(1j * (k* -LA(:,:,j+3))));
        DBAs = (DBAs + KBA1(:,:,j) * exp(1j * (k* FA(:,:,j))) + KBA3(:,:,j) * exp(1j * (k* TA(:,:,j))) +  ...
                KBA4(:,:,j) * exp(1j * (k* LA(:,:,j))) + KBA4(:,:,j+3) * exp(1j * (k* LA(:,:,j+3))));
    D(1:3, 4:6) = -DABs/sqrt(MN)/sqrt(MB);
    D(4:6, 1:3) = -DBAs/sqrt(MN)/sqrt(MB);
    D(1:3, 1:3) = (sum(KAB1,3)+sum(KAA,3)+sum(KAB3,3)+sum(KAB4f,3)+sum(KAB4s,3)-DAAs)/MB;
    D(4:6, 4:6) = (sum(KBA1,3) + sum(KBB,3) + sum(KBA3,3) + sum(KBA4f,3) + sum(KBA4s,3) - DBBs)/MN;
    D1(1:3, 4:6) = -DABs/sqrt(MN)/sqrt(MB1);
    D1(4:6, 1:3) = -DBAs/sqrt(MN)/sqrt(MB1);
    D1(1:3, 1:3) = (sum(KAB1,3)+sum(KAA,3)+sum(KAB3,3)+sum(KAB4f,3)+sum(KAB4s,3)-DAAs)/MB1;
    D1(4:6, 4:6) = (sum(KBA1,3) + sum(KBB,3) + sum(KBA3,3) + sum(KBA4f,3) + sum(KBA4s,3) - DBBs)/MN;
    D2(1:3, 4:6) = -DABs/sqrt(MN1)/sqrt(MB);
    D2(4:6, 1:3) = -DBAs/sqrt(MN1)/sqrt(MB);
    D2(1:3, 1:3) = (sum(KAB1,3)+sum(KAA,3)+sum(KAB3,3)+sum(KAB4f,3)+sum(KAB4s,3)-DAAs)/MB;
    D2(4:6, 4:6) = (sum(KBA1,3) + sum(KBB,3) + sum(KBA3,3) + sum(KBA4f,3) + sum(KBA4s,3) - DBBs)/MN1;
    D3(1:3, 4:6) = -DABs/sqrt(MN1)/sqrt(MB1);
    D3(4:6, 1:3) = -DBAs/sqrt(MN1)/sqrt(MB1);
    D3(1:3, 1:3) = (sum(KAB1,3)+sum(KAA,3)+sum(KAB3,3)+sum(KAB4f,3)+sum(KAB4s,3)-DAAs)/MB1;
    D3(4:6, 4:6) = (sum(KBA1,3) + sum(KBB,3) + sum(KBA3,3) + sum(KBA4f,3) + sum(KBA4s,3) - DBBs)/MN1;
    end
    [t,w]=eig(D);
    wg=sort(diag(w));
    [t,w]=eig(D1);
    wg1=sort(diag(w));
    [t,w]=eig(D2);
    wg2=sort(diag(w));
    [t,w]=eig(D3);
    wg3=sort(diag(w));
    result(i,:)=[i;sqrt(wg)/(299792458*100)/2/pi];
    result1(i,:)=[i;sqrt(wg1)/(299792458*100)/2/pi];
    result2(i,:)=[i;sqrt(wg2)/(299792458*100)/2/pi];
    result3(i,:)=[i;sqrt(wg3)/(299792458*100)/2/pi];
end


xlswrite('C:\Users\xz\Desktop\phonon\phonon.xlsx',result,'phonon0')
xlswrite('C:\Users\xz\Desktop\phonon\phonon.xlsx',result1,'phonon1')
xlswrite('C:\Users\xz\Desktop\phonon\phonon.xlsx',result2,'phonon2')
xlswrite('C:\Users\xz\Desktop\phonon\phonon.xlsx',result3,'phonon3')