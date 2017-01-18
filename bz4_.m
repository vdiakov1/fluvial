BZ_lib % functions are defined there


tf = 400; dt = 0.01; T= 0:dt:tf ; %tf = 400; dt = 0.004; T= 0:dt:tf ; % The time span
yo = [4 0.6 4]'; %[4 0.5 4]'; %[4 1.1 4]'; % Initial conitions
yss = steadystate(); J0 = joreg(yss,0); [A,B]=eig(J0) %eigvals = eig(J0); 
display("Steady state and Jacobian eigen values:")
[yss,diag(B),yporeg(yss,0)]
% solve ODEs:
[y, ISTATE, MSG] = lsode('yporeg',yo,T); % matlab: [t y] = ode15s('yporeg',[0 tf],yo);

[a,b,c] = callabc(); titlu = ['a=',num2str(a),', b=',num2str(b), ', c=',num2str(c)] ;
 plot_it(T, y, titlu);
  figure(2); hold on
   plot3(log(y(:,1)),log(y(:,2)),log(y(:,3)),'g--')
    xlabel('log(y1)');ylabel('log(y2)');zlabel('log(y3)')
     title(titlu)
Ly=floor(length(y)/2);y2max=max(y(Ly:(Ly*2),2)); imax = find(y(Ly:(Ly*2),2)==y2max)+Ly;
ystart = y(imax,:)';
yc1=fsolve(@gradconfluence,ystart)
figure(2); hold on
ys=[yc1';yc1']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'o')
ys=[yss';yss']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'x')

% build zero-lateral-gradient-of-confluence curve
% larger steps:
yc=zerogradfinder(ystart);[congrad(yc)'; yc']
z=yc'; scale=1;
while (size(z,1)<16)%&(validzerograd(z(size(z,1),:)'));
p=normalize(vectorproduct(yc-yss,A(:,1))); p' ; while (min(yc)>1e-5)*(min(p)<-1e-12)*(min(yc+p)<0) ; p = p/2.0 ;endwhile
yc=zerogradfinder(yc+p);[congrad(yc)'; yc'] ;
%if validzerograd(yc) ;;else;scale /= 2.0; yc -= scale*p; yc=zerogradfinder(yc); endif
if validzerograd(yc); z=[z;yc'];
else;display(size(z,1));display([congrad(yc)'; yc']);scale /= 2.0; yc -= scale*p; yc=zerogradfinder(yc) ;endif ; %p = normalize(yc - z(size(z,1)-1,:)')
endwhile
plot3(log(z(:,1)),log(z(:,2)),log(z(:,3)),'b.')

v=0*z(:,1);for i=1:size(z,1);v(i)=abscongrad(z(i,:)');endfor % imprecision of congrad == 0
find(v>0.1)'

zs=[0 0 0];for fr=[0.2 0.15, -0.1, -0.2];yc=zerograd((1-fr)*yc1+fr*yss);display([[fr 0 cosangle(yc-yc1,yss-yc)];congrad(yc)'; yc']);zs=[yc';zs];endfor;zs=zs(1:4,:);
fr=0.1;yc=zerogradfinder((1-fr)*yc1+fr*yss);display([[fr 0 cosangle(yc-yc1,yss-yc)];congrad(yc)'; yc']);
yc3 = yc; % keep this direction % 1.46518    1.94609   13.15287
frzz=[-2.0 -1.5 -1.25 -1.0 0.75 1 1.25 1.5 2];
zz=[0 0 0];for fr=frzz;yc=zerograd((1-fr)*yc1+fr*yc3);zz=[yc';zz];endfor;zz=zz(1:size(zz,1)-1,:);
zx=z(1:15,:);
z15=z(15,:)'; z16=z(16,:)';
frzx=[0.1,0.2,0.4,0.6];
for fr=frzx;yc=zerograd((1-fr)*z15+fr*z16);zx=[yc';zx];endfor

figure; hold on
   plot3(log(y(:,1)),log(y(:,2)),log(y(:,3)),'g--')
    xlabel('log(y1)');ylabel('log(y2)');zlabel('log(y3)')
     title(titlu)
ys=[yc1';yc1']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'o')
ys=[yss';yss']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'x')
plot3(log(zx(:,1)),log(zx(:,2)),log(zx(:,3)),'b.')
plot3(log(zs(:,1)),log(zs(:,2)),log(zs(:,3)),'r.')
plot3(log(zz(:,1)),log(zz(:,2)),log(zz(:,3)),'k.')

A1=[0.017757   0.147189  -0.017607;
   0.068990  -0.199965   0.091805;
  -0.997459  -0.968684   0.995621];

% solve ODEs:
y0=[2.7219    2.3557   14.0128]';%=zx(14,:)'+0.4*A1(:,1);
[y1, ISTATE, MSG] = lsode('yporeg',y0,0:dt/10:5);
y0=[2.6971    2.2591   15.8092]'%=zx(14,:)'-1*A1(:,1); % two 'joining' solutions in a couloir
[y2, ISTATE, MSG] = lsode('yporeg',y0,0:dt:3); % this run is just to set up y10 and y11
y2(size(y2,1),:) %    1.6886    2.3703   13.6779
[y2, ISTATE, MSG] = lsode('yporeg',y0,0:dt/10:3);
y0=[1.2    2.26   15.8]';%=zx(14,:)'+0.4*A1(:,1);
[y3, ISTATE, MSG] = lsode(@yporeg,y0,0:dt/10:5);
y10=zerograd(y2(size(y2,1),:)'); y10' %   1.6881    2.3710   13.6779 
y11=zerograd(y2(150,:)'); y11' % 1.7284    2.2912   15.1763
fr=-0.3975;yc=zerograd((1-fr)*y11+fr*y10);display([[fr 0 cosangle(yc-y11,y10-yc)];congrad(yc)'; yc']);
y12 = yc; % 1.74604    2.25920   15.77197
fr=1.1;yc=zerograd((1-fr)*y11+fr*y10);display([[fr 0 cosangle(yc-y11,y10-yc)];congrad(yc)'; yc']);
y01=yc; %  1.68620    2.37501   13.52782
zy=[y12';y11';y10';y01'];
[y2, ISTATE, MSG] = lsode('yporeg',y0,0:dt/10:5);
yc=zerograd(y2(4001,:)');display([[fr 0 cosangle(yc-y11,y10-yc)];congrad(yc)'; yc']);
zy(5,:)=yc';% 1.66575    2.41905   12.67933
yc=zerograd(y2(5001,:)');display([[fr 0 cosangle(yc-y11,y10-yc)];congrad(yc)'; yc']);
zy(6,:)=yc';% 1.64687    2.46228   11.68145
%{yc=next0grad(4.0,zy);%s=size(zy,1);z1=zy(s,:)';z2=zy(s-1,:)';fr=4.0;yc=zerograd((1-fr)*z2+fr*z1);display([[fr 0 cosangle(yc-z1,yc-z2)];congrad(yc)'; yc']);
zy=[zy;yc']; % 1.68066    2.38654   13.07749
yc=next0grad(3.305,zy);
zy=[zy;yc']; % 1.66887    2.41205   12.03944
%}
figure; hold on % this figure with calcs from bzlib (not bzlib4)
plot3(log(y1(:,1)),log(y1(:,2)),log(y1(:,3)),'g--')
    xlabel('log(y_1)');ylabel('log(y_2)');zlabel('log(y_3)')
     title(titlu)
plot3(log(y2(:,1)),log(y2(:,2)),log(y2(:,3)),'g--')
plot3(log(y3(:,1)),log(y3(:,2)),log(y3(:,3)),'g--')
plot3(log(zy(:,1)),log(zy(:,2)),log(zy(:,3)),'r.')

 jac = lateraljacobian(@congrad,zy(1,:)')
[A,B] = eig(jac)
cosangle(y2(1,:)'-zy(1,:)',A(:,1)) % 0.9 => mostly from the strongest eig.vector direction
cosangle(y2(5001,:)'-zy(1,:)',A(:,1)) % 0.01
congrad(zy(1,:)')' %  -11.1152   24.9403    1.5215 => much smaller than in the vicinity (see below):
congrad(zy(1,:)'+0.01*(y2(1,:)-zy(1,:))')' %  -1449.13   -981.48    590.57
congrad(zy(1,:)'+0.001*(y2(1,:)-zy(1,:))')'%  -272.7924  -123.0571     8.6476

zyd=[];for i=1:size(zy,1); zyd = [zyd depth(zy(i,:)')]; endfor; zyd
zxd=[];for i=1:size(zx,1); zxd = [zxd depth(zx(i,:)')]; endfor; zxd
zsd=[];for i=1:size(zs,1); zsd = [zsd depth(zs(i,:)')]; endfor; zsd
zzd=[];for i=1:size(zz,1); zzd = [zzd depth(zz(i,:)')]; endfor; zzd

[y3, ISTATE, MSG] = lsode('yporeg',y0,0:dt:10);

fr=2.;yc=zerograd((1-fr)*zx(2,:)'+fr*zx(1,:)');display([[fr depth(yc) cosangle(yc-zx(2,:)',zx(1,:)'-yc)];congrad(yc)';yc']);
zm=yc';% 2.71476    0.56848   15.13566
fr=2.;yc=zerograd((1-fr)*zx(1,:)'+fr*zm(1,:)');display([[fr depth(yc) cosangle(yc-zx(1,:)',zm(1,:)'-yc)];congrad(yc)';yc']);
zm=[zm;yc'];% 2.7150e+000  3.7006e-001  1.5153e+001
fr=2.;s=size(zm,1);ym=zm(s,:)';yn=zm(s-1,:)';yc=zerograd((1-fr)*yn+fr*ym);display([[fr depth(yc) cosangle(yc-yn,ym-yc)];congrad(yc)';yc']);
zm=[zm;yc'];% 2.71523    0.17139   15.17000
fr=2.;s=size(zm,1);ym=zm(s,:)';yn=zm(s-1,:)';yc=zerograd((1-fr)*yn+fr*ym);display([[fr depth(yc) cosangle(yc-yn,ym-yc)];congrad(yc)';yc']);
zm=[zm;yc'];% 1.8384e+001  1.1530e-001  1.1260e+001
fr=2.;s=size(zm,1);ym=zm(s,:)';yn=zm(s-1,:)';yc=zerograd((1-fr)*yn+fr*ym);display([[fr depth(yc) cosangle(yc-yn,ym-yc)];congrad(yc)';yc']);
zm=[zm;yc'];% 3.4054e+001  9.9999e-002  7.3521e+000
fr=2.5;s=size(zm,1);ym=zm(s,:)';yn=zm(s-1,:)';yc=zerograd((1-fr)*yn+fr*ym);display([[fr depth(yc) cosangle(yc-yn,ym-yc)];congrad(yc)';yc']);
zm=[zm;yc'];% 
% more details needeed between zm(2,:) and zm(3,:); also zm(3,:)-zm(4,:)
yc=next0grad(0.95,zm(1:3,:));% 2.7152e+000  1.8136e-001  1.5169e+001
zm=[zm(1:2,:);yc';zm(3:size(zm,1),:)];
yc=next0grad(2,zm(1:4,:));% 2.7152e+000  1.6138e-001  1.5171e+001
zm=[zm(1:4,:);yc';zm(5:size(zm,1),:)];
yc=next0grad(2,zm(1:5,:));% 2.7152e+000  1.5119e-001  1.5172e+001
zm=[zm(1:5,:);yc';zm(6:size(zm,1),:)];
yc=next0grad(3,zm(1:6,:));% 2.7153e+000  1.3128e-001  1.5174e+001
zm=[zm(1:6,:);yc';zm(7:size(zm,1),:)];
yc=next0grad(2,zm(1:7,:));% 2.71530    0.11176   15.17517
zm=[zm(1:7,:);yc';zm(8:size(zm,1),:)];
yc=next0grad(1.75,zm(1:8,:));% 2.7153e+000  1.0000e-001  1.5176e+001
zm=[zm(1:8,:);yc';zm(9:size(zm,1),:)];
%{yc=next0grad(2,zm(1:9,:));% 2.71549    0.10000   15.18724
zm=[zm(1:9,:);yc';zm(10:size(zm,1),:)];
yc=next0grad(2,zm(1:10,:));% 2.71568    0.10013   15.19798
zm=[zm(1:10,:);yc';zm(11:size(zm,1),:)];
yc=next0grad(3,zm(1:11,:));% 2.71605    0.10042   15.21944
zm=[zm(1:11,:);yc';zm(12:size(zm,1),:)];
yc=next0grad(3,zm(1:12,:));% 2.71679    0.10084   15.26245
zm=[zm(1:12,:);yc';zm(13:size(zm,1),:)];
yc=next0grad(4,zm(1:13,:));% 2.71901    0.10217   15.39152
zm=[zm(1:13,:);yc';zm(14:size(zm,1),:)];
yc=next0grad(4,zm(1:14,:));% 2.7257e+000  1.0612e-001  1.5779e+001
zm=[zm(1:14,:);yc';zm(15:size(zm,1),:)];
yc=next0grad(4,zm(1:15,:));% 2.7456e+000  1.1806e-001  1.6940e+001
zm=[zm(1:15,:);yc';zm(16:size(zm,1),:)];
yc=next0grad(4,zm(1:16,:));% 2.80553    0.15360   20.42600
zm=[zm(1:16,:);yc';zm(17:size(zm,1),:)];
yc=next0grad(2,zm(1:17,:));% 2.86542    0.18954   23.91106
zm=[zm(1:17,:);yc';zm(18:size(zm,1),:)];
yc=next0grad(2,zm(1:18,:));% 2.92529    0.22634   27.39533
zm=[zm(1:18,:);yc';zm(19:size(zm,1),:)];
%}
% continue from zm(9,:) to a different direction (towards zm(size(zm,1)-1,:) )
s=size(zm,1);z1=zm(s,:)';z2=zm(s-1,:)';zm=zm(1:9,:); % z2'=[18.38447    0.11530   11.26012] 
yc=next0grad(0.05,[zm;z2']);% 3.4988e+000  1.0054e-001  1.4980e+001
zm=[zm;yc'];
yc=next0grad(0.25,[zm;z2']);% 7.2202e+000  1.0414e-001  1.4050e+001
zm=[zm;yc'];
yc=next0grad(0.33,[zm;z2']);% 1.0904e+001  1.1013e-001  1.3128e+001
zm=[zm;yc'];
yc=next0grad(0.5,[zm;z2']);% 1.4644e+001  1.1404e-001  1.2194e+001
zm=[zm;yc'];
yc=next0grad(0.65,[zm;z2']);% 1.7075e+001  1.1365e-001  1.1588e+001
zm=[zm;yc'];
yc=next0grad(0.5,[zm;z2']);% 1.7730e+001  1.1321e-001  1.1425e+001
zm=[zm;yc'];
zm=[zm;z2']; % end interpolating towards z2, start extrapolating; total 16 points in zm so far
%{ below is a different branch
yc=next0grad(3,zm);% 1.9693e+001  1.2226e-001  1.0925e+001
zm=[zm;yc'];
yc=next0grad(2,zm);% 2.1003e+001  1.4219e-001  1.0590e+001
zm=[zm;yc'];
yc=next0grad(2,zm);% 2.2312e+001  1.6204e-001  1.0258e+001
zm=[zm;yc']; % 19 points in zm
yc=next0grad(3,zm);% 24.92955    0.20736    9.59610
zm=[zm;yc']; 
yc=next0grad(2,zm);% 2.7548e+001  2.5375e-001  8.9296e+000
zm=[zm;yc']; 
yc=next0grad(1.5,zm);% 2.8857e+001  2.7032e-001  8.6021e+000
zm=[zm;yc']; 
yc=next0grad(3,zm);% 31.47439    0.28415    7.95307
zm=[zm;yc']; % 23
yc=next0grad(1.75,zm);% 33.43851    0.28968    7.44481
zm=[zm;yc'];
yc=next0grad(3,zm);% 37.36718    0.31042    6.43852
zm=[zm;yc']; % 25
yc=next0grad(1.67,zm);% 39.99939    0.33756    5.76424
zm=[zm;yc']; % 26
%}
yc2=next0grad(0.25,[zm(1:16,:);z1']);% 2.2302e+001  1.1003e-001  1.0283e+001
%yc3=next0grad(0.45,[zm(1:16,:);z1']);% 2.5436e+001  1.2173e-001  9.5024e+000
yc4=next0grad(0.65,[zm(1:16,:);z1']);% 2.8569e+001  1.0145e-001  8.7199e+000
yc5=next0grad(0.75,[zm(1:16,:);z1']);% 3.0136e+001  1.0649e-001  8.3318e+000
zm=[zm(1:16,:);yc2';yc4';yc5';z1'];%yc3';
yc=next0grad(2,zm);% 3.7971e+001  1.0000e-001  6.3708e+000
zm=[zm;yc']; % 22
yc=next0grad(1.1,zm);% 3.8363e+001  1.0000e-001  6.2727e+000
zm=[zm;yc']; % 23
yc=next0grad(4.9,zm);% 3.9890e+001  1.1122e-001  5.8938e+000
zm=[zm;yc']; % 24
z4=[2.65 0.1 5.51]';% approx. expected coords of the 2nd 'corner'
yc2=next0grad(0.15,[zm;z4']);% 3.4304e+001  1.0852e-001  5.8383e+000
yc3=next0grad(0.23,[zm;z4']);% 3.1325e+001  1.0494e-001  5.8074e+000
yc5=next0grad(0.525,[zm;z4']);% 2.0339e+001  1.0756e-001  5.6916e+000
yc6=next0grad(0.625,[zm;z4']);% 1.6615e+001  1.0597e-001  5.6539e+000
yc7=next0grad(0.755,[zm;z4']);% 1.1774e+001  1.0153e-001  5.6065e+000
yc8=next0grad(0.845,[zm;z4']);% 8.42226   0.10376   5.56980
yc91=next0grad(1,[zm;z4']);% 2.6500e+000  1.0000e-001  5.5101e+000
zm=[zm;yc2';yc3';yc5';yc6';yc7';yc8';yc91'];
% new direction: towards z(1,:):
yc2=next0grad(0.2,[zm;z(1,:)]);% 2.6496e+000  3.8290e-001  5.5177e+000
yc4=next0grad(0.44,[zm;z(1,:)]);% 2.6491e+000  7.1710e-001  5.5280e+000
yc6=next0grad(0.64,[zm;z(1,:)]);% 2.6487e+000  1.0089e+000  5.5304e+000
yc8=next0grad(0.84,[zm;z(1,:)]);% 2.6483e+000  1.3039e+000  5.5389e+000
zm=[zm;yc2';yc4';yc6';yc8'];
plot3(log(zm(:,1)),log(zm(:,2)),log(zm(:,3)),'b.')
plot3(log(zm([18 19 21 23 24 25 27],1)),log(zm([18 19 21 23 24 25 27],2)),log(zm([18 19 21 23 24 25 27],3)),'w.')% these points are off
zn=[zm(1:17,:);zm(20,:);zm(22,:);zm(26,:);zm(28:36,:)];% exclude the 'off' points
plot3(log(zn(:,1)),log(zn(:,2)),log(zn(:,3)),'b.')
plot3(log(zx(:,1)),log(zx(:,2)),log(zx(:,3)),'b.')


figure(1); hold on
%ys=[yc1';yc1']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'wo')
plot3(log(zn(:,1)),log(zn(:,2)),log(zn(:,3)),'b.') % show the missing part of the 'boot'
%{
>> [[1:size(zn,1)]',zn]
    1.00000    2.71476    0.56848   15.13566
    2.00000    2.71499    0.37006   15.15284
    3.00000    2.71521    0.18136   15.16909
    4.00000    2.71523    0.17139   15.17000
    5.00000    2.71524    0.16138   15.17089
    6.00000    2.71525    0.15119   15.17173
    7.00000    2.71527    0.13128   15.17352
    8.00000    2.71530    0.11176   15.17517
    9.00000    2.71531    0.10000   15.17649
   10.00000    3.49877    0.10054   14.98045
   11.00000    7.22019    0.10414   14.05031
   12.00000   10.90440    0.11013   13.12805
   13.00000   14.64444    0.11404   12.19363
   14.00000   17.07546    0.11365   11.58784
   15.00000   17.72996    0.11321   11.42501
   16.00000   18.38447    0.11530   11.26012
   17.00000   22.30176    0.11003   10.28290
   18.00000   30.13632    0.10649    8.33185
   19.00000   37.97094    0.10000    6.37085
   20.00000   31.32507    0.10494    5.80740
   21.00000   20.33917    0.10756    5.69160
   22.00000   16.61514    0.10597    5.65393
   23.00000   11.77388    0.10153    5.60653
   24.00000    8.42226    0.10376    5.56980
   25.00000    2.65000    0.10000    5.51012
   26.00000    2.64958    0.38290    5.51774
   27.00000    2.64909    0.71710    5.52804
   28.00000    2.64866    1.00895    5.53040
   29.00000    2.64829    1.30394    5.53895%}
%{
yc=next0grad(2,zx);% 
zp=yc';% 2.71287    0.39976   15.32494
yc=next0grad(1.2225,[zx;zp]);% 2.7127e+000  1.8528e-001  1.5382e+001
zp=[zp;yc'];
yc=next0grad(2,zp);% 2.7137e+000  9.9983e-002  1.5383e+001
zp=[zp;yc']; size(zp,1)

>
>> p=normalize(vectorproduct(z(size(z,1),:)'-yss,A(:,1))); p'
ans =   0.0015649  -0.9996078  -0.0279623
>> yc=zerogradfinder(z(size(z,1),:)'+p);congrad(yc)'
ans =  -5.6753e-005  1.2864e-002  -1.2784e-004
>> z=[z;yc'] ; z(size(z,1),:)
ans =    6.65988    0.10000   35.66479

z =[2.647878    1.521757    5.538916
    2.648454    2.507155    5.692412
    2.650744    3.430536    6.076129
    2.654578    4.239089    6.663407
    2.659713    4.894988    7.418422
    2.665861    5.367938    8.299333
    2.672705    5.638215    9.262538
    2.679920    5.701459   10.260642
    2.687185    5.557568   11.250326
    2.694195    5.216257   12.190872
    2.700675    4.706510   13.044844
    2.706396    4.032421   13.782182
    2.711157    3.229658   14.377534
    2.714819    2.328113   14.811741
    2.713844    1.364040   15.068330
    2.714978    0.367776   15.151661
    2.714938    0.118808   15.127883
    2.768622    0.100012   15.106813
    2.852987    0.099995   15.113823
    2.961857    0.099999   15.220030
    3.072051    0.099990   15.302763
    3.074136    0.100007   15.276478
    3.170028    0.099983   15.265804
    3.339981    0.099988   14.899221
    3.433573    0.100007   14.835373
    3.311143    0.099999   13.881715
    3.424021    0.099988   13.918764
    3.423927    0.099998   13.907468
    3.453004    0.100006   14.102931
    3.453797    0.099997   14.079717
    3.455653    0.100006   14.126653
    3.554049    0.100007   14.080781
    3.557313    0.099996   14.120409];

%fr=0.2;yc=zerogradfinder((1-fr)*yc1+fr*yss);display([[fr 0 cosangle(yc-yc1,yss-yc)];congrad(yc)'; yc']);
% for fr=[0.2 0.15, -0.1, -0.2], do the line above and store results in zs:
%zs=[1.25726    1.65343   15.52499; 1.94523    1.59395   15.02478; 3.65179    1.43716   13.73647; 3.98436    1.40847   13.51648];
zs=[0 0 0];for fr=[0.2 0.15, -0.1, -0.2];yc=zerograd((1-fr)*yc1+fr*yss);display([[fr 0 cosangle(yc-yc1,yss-yc)];congrad(yc)'; yc']);zs=[yc';zs];endfor;zs=zs(1:4,:);
   % this is the line that crosses z(:,:) and goes thru the max-confluence yc1 and the steadystate yss points
% now, find the 3rd line that has to go thru the intersection of the other two (z and zs) lines:
fr=0.1;yc=zerogradfinder((1-fr)*yc1+fr*yss);display([[fr 0 cosangle(yc-yc1,yss-yc)];congrad(yc)'; yc']);
yc3 = yc; % keep this direction % 1.46518    1.94609   13.15287
fr=2;yc=zerogradfinder((1-fr)*yc1+fr*yc3);display([[fr 0 cosangle(yc-yc1,yc3-yc)];congrad(yc)'; yc']);
%for fr=[-2.0 -1.5 -1.25 -1.0 0.75 1 1.25 1.5 2]
%frzz=[-2.0 -1.5 -1.25 -1.0 0.75 1 1.25 1.5 2];
%zz=[0 0 0];for fr=frzz;yc=zerograd((1-fr)*yc1+fr*yc3);zz=[yc';zz];endfor;zz=zz(1:size(zz,1)-1,:);
zz=[4.96903 0.69017 17.26256; 4.38029 0.90136 16.59006; 4.10066 1.00602 16.22648; 3.80980 1.10982 15.89371; 
    1.76295 1.83586 13.4911; 1.46518 1.9461 13.15287; 1.17318 2.05035 12.81; 0.88121 2.15497 12.46785; 0.29723 2.36384 11.78293];
% adjust the elements of z:
zx=z(1:15,:);
z15=z(15,:)'; z16=z(16,:)';
%fr=0.1;yc=zerogradfinder((1-fr)*z15+fr*z16);display([[fr 0 cosangle(yc-z15,z16-yc)];congrad(yc)'; yc']);
%zx=[zx;yc']; % and so forth with fr=[0.1,0.2,0.4,0.6]
frzx=[0.1,0.2,0.4,0.6];
for fr=frzx;yc=zerograd((1-fr)*z15+fr*z16);zx=[yc';zx];endfor
figure; hold on
   plot3(log(y(:,1)),log(y(:,2)),log(y(:,3)),'g--')
    xlabel('log(y1)');ylabel('log(y2)');zlabel('log(y3)')
     title(titlu)
ys=[yc1';yc1']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'o')
ys=[yss';yss']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'x')
plot3(log(zx(:,1)),log(zx(:,2)),log(zx(:,3)),'b.')
plot3(log(zs(:,1)),log(zs(:,2)),log(zs(:,3)),'r.')
plot3(log(zz(:,1)),log(zz(:,2)),log(zz(:,3)),'k.')

A1=[0.017757   0.147189  -0.017607;
   0.068990  -0.199965   0.091805;
  -0.997459  -0.968684   0.995621];
a3=A1(:,3); zx0=zx(15,:)'-zx(14,:)'; zz0=yc3-yc1; zs0 = zs(2,:)'-zs(3,:)'; % directions: the weakest eig.v. and 3 curves
display([cosangle(a3,zx0) cosangle(a3,zz0) cosangle(a3,zs0)]) % zs is most aligned with the weakest eig.v.

% solve ODEs:
y0=[2.7219    2.3557   14.4128]';%=zx(14,:)'+0.4*A1(:,1);
[y1, ISTATE, MSG] = lsode('yporeg',y0,0:dt/10:5);
y0=[2.6971    2.2591   15.8092]'%=zx(14,:)'-1*A1(:,1); % two 'joining' solutions in a couloir
[y2, ISTATE, MSG] = lsode('yporeg',y0,0:dt:3); % this run is just to set up y10 and y11
y2(size(y2,1),:) %    1.6886    2.3703   13.6779
y10=zerograd(y2(size(y2,1),:)'); y10' %   1.6888    2.3693   13.6778 
y11=zerograd(y2(150,:)'); y11' % 1.7268    2.2945   15.1855
[y2, ISTATE, MSG] = lsode('yporeg',y0,0:dt/10:3);
fr=-0.4;yc=zerograd((1-fr)*y11+fr*y10);display([[fr 0 cosangle(yc-y11,y10-yc)];congrad(yc)'; yc']);
y12 = yc; % 1.74413    2.26254   15.78852
fr=1.1;yc=zerograd((1-fr)*y11+fr*y10);display([[fr 0 cosangle(yc-y11,y10-yc)];congrad(yc)'; yc']);
y01=yc; %  1.68393    2.37981   13.52727
zy=[y12';y11';y10';y01'];
fr=1.6;yc=zerograd((1-fr)*y11+fr*y10);display([[fr 0 cosangle(yc-y11,y10-yc)];congrad(yc)'; yc']);
zy=[zy;yc']; % 1.66735     2.41581    12.77370
fr=2.7;yc=zerograd((1-fr)*y11+fr*y10);display([[fr 0 cosangle(yc-y11,y10-yc)];congrad(yc)'; yc']);
zy=[zy;yc']; % 1.60738     2.50544    11.10958
figure; hold on
   plot3(log(y(:,1)),log(y(:,2)),log(y(:,3)),'g--')
    xlabel('log(y1)');ylabel('log(y2)');zlabel('log(y3)')
     title(titlu)
ys=[yc1';yc1']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'o')
ys=[yss';yss']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'x')
plot3(log(y1(:,1)),log(y1(:,2)),log(y1(:,3)),'g--')
plot3(log(y2(:,1)),log(y2(:,2)),log(y2(:,3)),'g--')
plot3(log(zy(:,1)),log(zy(:,2)),log(zy(:,3)),'r.')

%}
