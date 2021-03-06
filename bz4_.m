BZ_lib % functions are defined there

% Oregonator solution and a couple of plots:
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

% examine fluvii:
z =[2.647878    1.521757    5.538916;    2.648454    2.507155    5.692412;
    2.650744    3.430536    6.076129;    2.654578    4.239089    6.663407;
    2.659713    4.894988    7.418422;    2.665861    5.367938    8.299333;
    2.672705    5.638215    9.262538;    2.679920    5.701459   10.260642;
    2.687185    5.557568   11.250326;    2.694195    5.216257   12.190872;
    2.700675    4.706510   13.044844;    2.706396    4.032421   13.782182;
    2.711157    3.229658   14.377534;    2.714819    2.328113   14.811741;
    2.713844    1.364040   15.068330;    2.714978    0.367776   15.151661];
zx=z(1:15,:);
z15=z(15,:)'; z16=z(16,:)';
frzx=[0.1,0.2,0.4,0.6];
for fr=frzx;yc=zerograd((1-fr)*z15+fr*z16);zx=[yc';zx];endfor
% how close are these points to fluvii:
zxcheck = from_fluvium(zx);

zz=[4.96903 0.69017 17.26256; 4.38029 0.90136 16.59006; 4.10066 1.00602 16.22648; 3.80980 1.10982 15.89371; 
    1.76295 1.83586 13.4911; 1.46518 1.9461 13.15287; 1.17318 2.05035 12.81; 0.88121 2.15497 12.46785];% 0.29723 2.36384 11.78293];
zzcheck = from_fluvium(zz);

zs=[1.25726    1.65343   15.52499; 1.94523    1.59395   15.02478; 3.65179    1.43716   13.73647; 3.98436    1.40847   13.51648];
zscheck = from_fluvium(zs);

znn = [
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
   29.00000    2.64829    1.30394    5.53895];
zy=znn(:,2:size(znn,2));
zycheck = from_fluvium(zy);

figure; hold on
   plot3(log(y(:,1)),log(y(:,2)),log(y(:,3)),'g--')
    xlabel('log(y1)');ylabel('log(y2)');zlabel('log(y3)')
     title(titlu)
ys=[yc1';yc1']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'o')
ys=[yss';yss']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'x')
plot3(log(zx(:,1)),log(zx(:,2)),log(zx(:,3)),'b.')
plot3(log(zy(:,1)),log(zy(:,2)),log(zy(:,3)),'b.')
plot3(log(zs(:,1)),log(zs(:,2)),log(zs(:,3)),'r.')
plot3(log(zz(:,1)),log(zz(:,2)),log(zz(:,3)),'k.')
