BZ_lib % load function definitions
% re-define parameter values
function [a,b,c] = callabc(); a=77.; b=0.161; c=1e-5; end % c=2e-2; end %

% Define time span:
tf = 400; dt = 0.004; T= 0:dt:tf ;
%tf = 120; dt = 0.01; T= 0:dt:tf ;

% Initial conitions:
yo = [4 0.6 4]'; %[4 0.5 4]'; %[4 1.1 4]'; %
yss = steadystate(); J0 = joreg(yss,0); eigvals = eig(J0); [A,B]=eig(J0)
display("Steady state and Jacobian eigen values:")
[yss,eigvals,yporeg(yss,0)]
% solve ODEs:
[y, ISTATE, MSG] = lsode('yporeg',yo,T); % matlab [t y] = ode15s('yporeg',[0 tf],yo);
[a,b,c] = callabc(); titlu = ['a=',num2str(a),', b=',num2str(b), ', c=',num2str(c)] ;
 plot_it(T, y, titlu);
  figure(2); hold on
   plot3(log(y(:,1)),log(y(:,2)),log(y(:,3)),'g--')
    xlabel('log(y_1)');ylabel('log(y_2)');zlabel('log(y_3)')
     title(titlu)

     Ly=floor(length(y)/2);y2max=max(y(Ly:(Ly*2),2)); imax = find(y(Ly:(Ly*2),2)==y2max)+Ly;
ystart = y(imax,:)';
yc1=fsolve(@gradconfluence,ystart); yc1'
figure(2); hold on
ys=[yc1';yc1']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'o')
ys=[yss';yss']; plot3(log(ys(:,1)),log(ys(:,2)),log(ys(:,3)),'x')

y1s = exp([0:0.1:12]); y2s=exp([-6:0.1:8]); y3s=exp([-2:0.1:12]);
n1=length(y1s); n2=length(y2s); n3=length(y3s);
cs = zeros(n1,n2,n3); cx=[]; cy=cs; % initialize confluence values grid
% calculate confluence on the y1s x y2s x y3s grid:
for i=1:n1; for j=1:n2; for k=1:n3; cs(i,j,k)=confluence([y1s(i),y2s(j),y3s(k)]') ;endfor;endfor;endfor

cp=max1D(y1s,y2s,y3s,cs);% calculate the one-directional 'sticking planes' on the grid
yc=entrained(y); % trajectory points that are entrained in fluvii

  figure; hold on
   plot3(log(cp(:,1)),log(cp(:,2)),log(cp(:,3)),'g.', 'markersize',1)
    plot3(log(y(:,1)),log(y(:,2)),log(y(:,3)),'r-','linewidth',0.99)
     xlabel('log(y_1)');ylabel('log(y_2)');zlabel('log(y_3)')
      title(titlu)
  plot3(log(yc(:,1)),log(yc(:,2)),log(yc(:,3)),'rx', 'markersize',7)

y1t = exp([11.2:0.005:11.6]); y2t=exp([-5:0.05:-2]); y3t=exp([7:0.05:10]);
 n1=length(y1t); n2=length(y2t); n3=length(y3t);
  ct = zeros(n1,n2,n3);
   for i=1:n1; for j=1:n2; for k=1:n3; ct(i,j,k)=confluence([y1t(i),y2t(j),y3t(k)]') ;endfor;endfor;endfor
    cpt=max1D(y1t,y2t,y3t,ct);% 'sticking planes' in more detail
figure; hold on
 plot3(log(cpt(:,1)),log(cpt(:,2)),log(cpt(:,3)),'g.', 'markersize',2)
  xlabel('log(y_1)');ylabel('log(y_2)');zlabel('log(y_3)')
   axis([11.2,11.6,-5,-2,7,10])
    plot3(log(y(:,1)),log(y(:,2)),log(y(:,3)),'r.','markersize',3)
     title(titlu)
