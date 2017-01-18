% Oregonator ODEs (following http://www4.ncsu.edu/eos/users/w/white/www/white/ma302/less15.PDF):
%{ Brief description of functions below:
callabc() = parameter values for use in yporeg()
box_it() - not used; attempts to keep all operations within a box
fx(y) = the 'y1' (or 'x') component of derivatives yporeg(), i.e. dy/dt
        y is a column vector with 3 components (here and throughout this module)
fy() and fz() - same as fx(), only for 'y2' (or 'y') and 'y3' (or 'z'), respectively
fsquare(y) = the square of yporeg(y)
yporeg(y) = the vector of derivatives dy/dt as defined by the Oregonator model
steadystate() = the stationary point (where dy/dt = 0)
deriv(fn,y,dy) = the 'derivative' of fn at point y for increment dy (y & dy are column vectors)
ascent(fn,y) = the gradient of scalar function fn
divrgnce(fn,y) = the divergence of fn at point y
joreg(y,t) = closed-form-calculated Jacobian of yporeg() at point y (t is not used)
jacobian(fn,y) = the Jacobian of fn at point y
lateraljacobian(fn,y) = Jacobian of vector-function fn, with removed component along yporeg()
mid_eig(y) = medium eigen value for the Jacobian of congrad()
max_eig(y) = maximum eigen value for the Jacobian of congrad()
confluence(y) = negative lateral (i.e. yporeg-component removed) divergence of yporeg()
mconfluence(y) = - confluence(y)
gradconfluence(y) = gradient of confluence(y)
congrad(y) = lateral gradient (again, 'lateral' means that the yporeg-component is removed) of confluence
gradgradconfluence - not used
trackzerocongrads(y) = a guesstimate for another fluvial point assuming y belongs to a fluvium
                       (fluvii consist of points that maximally attract lateral trajectories
                        and are defined by lateralgradient(confluence)=0  )
                       Perhaps, this function can be improved by using lateraljacobian instead of jacobian
buildzerocurve(yo,dy,n) - finds n fluvial points starting with yo+dy
                          perhaps, can be improved by replacing @gradconfluence with @congrad
zerogradfinder(y) = fluvial locus using y as the starting point
zerograd(y) = zerogradfinder(y) improved by replacing jacobian() with lateraljacobian()
entrained(y) - scans ODE solution y looking for points that are close to fluvii (i.e. are 'entrained')
direction(y,p) = the search direction for the fluvium (i.e. the zero-congrad curve)
abscongrad(y) = the absolute value of congrad(y)
nextzerograd(y,p) = next fluvial (or, more precisely, zero-congrad) point at distance p from y
thefirsttwo(y) = two fluvial points, uses y as starting point
zerogradcurve(y,N) - finds N points where congrad()=0
validzerograd(y) - checks if y qualifies as a lateral attractor
next0grad(fr,zy) = next point with congrad()=0, uses the last two points in zy to extrapolate
                   the starting point with scaling factor fr
depth(y) = lateral attraction strength estimate
dispersiveness(y) - not used
max1D(y1s,y2s,y3s,cs) - from cs values defined on the grid y1s x y2s x y3s,
                        finds the 'grid-local' maxima in either 'x', 'y' or 'z' direction
nvrs(j) = inverse of j excluding the smallest eigen-value vector
cosangle(y1,y2) = cosine of the angle between vectors y1 and y2
distance(y1,y2) = distance between points y1 and y2
vectorproduct(x,y) = 3-D vector product of x and y
normalize(x) = vector x normalized to length 1
%}
function [a,b,c] = callabc(); a=77.; b=0.161; c=2e-2; end % c=1e-5; end %
function res = box_it(x,r1,r2);res=10000*((r1-x)^2*(x<r1)+(x-r2)^2*(x>r2)) ;end % penalty for off-limits 
function y1=fx(y,a,c); y1 = a*(y(2)-y(1)*y(2)+y(1)-c*y(1)*y(1)); end%
function y2=fy(y,a)  ; y2 = (-y(2) -y(1)*y(2)+y(3))/a ; end%
function y3=fz(y,b)  ; y3 = b*(y(1) - y(3)) ; end% 
function f2=fsquare(y,a,b,c) ; f2 = fx(y,a,c)^2 + fy(y,a)^2 + fz(y,b)^2 ;end
function dy=yporeg(y,t); [a,b,c]=callabc(); dy= [fx(y,a,c) fy(y,a) fz(y,b)]'; end%/sqrt(fsquare(y,a,b,c))
function y=steadystate(); [a,b,c]=callabc(); y1=(-c+sqrt(8*c+c*c))/(2*c); y2=y1/(1+y1); y3=y1; y=[y1,y2,y3]'; end
% gradient, divergence:
function df = deriv(fn,y,dy); dx=1.0/sqrt(dy'*dy); df = (fn(y+0.5*dy)-fn(y-0.5*dy))*dx; end % derivative along a vector dy
function grad = ascent(fn,y);dy=diag(0*y+1e-6);grad=0*y'; for n=1:length(y); grad(n)=deriv(fn,y,dy(:,n)) ;endfor; end % calculates gradient for scalar functions of vectors
function d = divrgnce(fn,y); jac = jacobian(fn,y); d = trace(jac); end
% Jacobian:
function jac=joreg(y,t); [a,b,c] = callabc(); jac=[a*(1-y(2)-2*c*y(1)) a*(1-y(1)) 0; -y(2)/a -(1+y(1))/a 1/a; b 0 -b]; end
function jac = jacobian(fn,y);dy=diag(0*y+1e-6);jac=0*dy; for n=1:length(y); jac(:,n)=deriv(fn,y,dy(:,n)) ;endfor; end % calculates Jacobian for vector-functions of vectors
function lj = lateraljacobian(fn,y);j=jacobian(fn,y);f=yporeg(y);lj=j-(j*f)*f'/sqrt(f'*f) ;end
function md=mid_eig(y);j=lateraljacobian(@congrad,y);e=eig(j);m=min(abs(e));ma=max(abs(e));n=find(abs(e)>m);md=-max(e(n));end%md=-mean(e(n));end%
function md=max_eig(y);j=lateraljacobian(@congrad,y);e=real(eig(j));md=max(e);end
% confluence:
function cf = confluence(y);fn=@yporeg;f=fn(y,0);f2=f'*f; dfidyj=jacobian(fn,y); d=trace(dfidyj); cf=f'*dfidyj*f/f2-d; end% lateral divergence (negative)
function mcf = mconfluence(y); mcf = - confluence(y) ;end % negative confluence (lateral divergence)
% confluence gradients:
function grad = gradconfluence(y);fn=@confluence;grad = ascent(fn,y) ;end % gradient of confluence - a separate function for fsolve
function cg = congrad(y);vf=@yporeg;sf=@confluence;f=vf(y,0);grad=ascent(sf,y)';cg=grad-f*(grad'*f)/(f'*f) ;end% lateral gradient of confluence; vector f-n vf (e.g.yporeg), scalar f-n sf (e.g.confluence)
function grads = gradgradconfluence(y);grads=jacobian(@congrad,y)' ;end % gradients (columns) of lateral gradients of confluence
% track fluvial loci:
function yc = trackzerocongrads(yd); fn=@congrad; grads=jacobian(fn,yd); vp=-vectorproduct(grads(:,1),grads(:,2)); yc1=yd+vp*0.01/sqrt(vp'*vp); yc=fsolve(fn,yc1); end
function zcd = buildzerocurve(yo,dy,nmax);fn=@gradconfluence;zcd=zeros(nmax,length(yo));zcd(1,:)=fsolve(fn,yo)';zcd(2,:)=fsolve(fn,zcd(1,:)'+dy)';for n=3:nmax;zcd(n,:)=fsolve(fn,2*zcd(n-1,:)'-zcd(n-2,:)')';endfor ;end
function y0=zerogradfinder(y);fn=@congrad;cg=fn(y);j=jacobian(fn,y)';y0=y-nvrs(j')*cg ;end % linear approximation for finding zero-gradient of confluence 
function y0=zerograd(y);fn=@congrad;cg=fn(y);j=lateraljacobian(fn,y)';y0=y-nvrs(j')*cg ;end % =zerogradfinder() with lateraljacobian() replacing jacobian()
function yc=entrained(y);yc=[]; for i=10000:1000:size(y,1);y0=y(i,:)';y1=zerograd(y0);if distance(y0,y1)/distance(y0,0*y0)<1e-3;yc=[yc;y1'];endif;endfor; end
function d=direction(y,p);j=gradgradconfluence(y);d=vectorproduct(j(:,1),j(:,2)) ;if cosangle(d,p)<0;d=-d;endif;end % find the search direction for the zero-congrad curve
function x=abscongrad(y);cg = congrad(y); x=sqrt(sum(cg.*cg)/sum(y.*y)) ; end %
function y1=nextzerograd(y,p);d=normalize(direction(y,p));y1=zerogradfinder(y+0.1*d);y2=zerogradfinder(y1);if abscongrad(y2)<abscongrad(y1);y1=y2;endif ;end % find next zero-congrad point
function zcd=thefirsttwo(y);y0=zerogradfinder(y);y0=zerogradfinder(y0);y1=nextzerograd(y0,[0 0 1e-6]');zcd=[y0';y1'] ;end % build the first two rows of zcd (i.e. zero-confluence-gradient curve)
function zcd=zerogradcurve(y,N);zcd=thefirsttwo(y); while size(zcd,1)<N;m=size(zcd,1);y0=nextzerograd(zcd(m,:)',zcd(m,:)'-zcd(m-1,:)');zcd=[zcd; y0'];endwhile;end
function is = validzerograd(y);m=min(y);ab=abscongrad(y);mi=max(abs(imag(i)));if (m<0)|(ab>0.5)|(mi>1e-3);is=false;else;is=true;endif ;end
function yc = next0grad(fr,zy);s=size(zy,1);z1=zy(s,:)';z2=zy(s-1,:)';yc=zerograd((1-fr)*z2+fr*z1);display([[fr depth(yc) cosangle(yc-z1,yc-z2)];congrad(yc)'; yc']);end
function d = depth(y);c=confluence(y);m=mid_eig(y);if (c>0)&(m>0);d=c^1.5/m^0.5;else;d=0;endif;d=d/sqrt(sum(y.*y)) ;end%estimate couloir depth
function d = dispersiveness(y);mc=-confluence(y);m=max_eig(y);if (mc>0)&(m>0);d=c^1.5/m^0.5;else;d=0;endif;d=d/sqrt(sum(y.*y)) ;end%estimate dispersiveness (opposite of depth)
% calculate one-directional 'sticking planes' on the grid
function cp=max1D(y1s,y2s,y3s,cs, cutoff = 0); cp =[];n1=length(y1s)-1;n2=length(y2s)-1;n3=length(y3s)-1;
 for i=2:n1; for j=2:n2; for k=2:n3; css=cs(i,j,k); ys = [y1s(i),y2s(j),y3s(k),css];cxm = css-cs(i-1,j,k); cxp = css-cs(i+1,j,k);
  cym=css-cs(i,j-1,k); cyp=css-cs(i,j+1,k);czm = css-cs(i,j,k-1); czp = css-cs(i,j,k+1);
   condition = ((cxm+cxp>0)*(cxm*cxp>0)+(cym+cyp>0)*(cym*cyp>0)+(czm+czp>0)*(czm*czp>0)) ;
    if condition*(css >= cutoff); cp = [cp;ys]; endif %*( css > 0 )
     endfor;endfor;endfor;end
% miscellaneous f-ns:
function B1=nvrs(j);[A B]=eig(j);b=diag(B);ab=abs(b);m=min(ab);if m<0.00001*max(ab);n=find(ab==m);b(n)=1;C=diag(1.0./b);C(n,n)=0;else;C=inv(B);endif;B1=real(A*C*inv(A)) ;end % inverse of j excluding the smallest eigen value vector
function ca=cosangle(vec1,vec2);ca=vec1'*vec2/sqrt((vec1'*vec1)*(vec2'*vec2)) ;end
function d=distance(vec1,vec2);v=vec1.-vec2; d=sqrt(sum(v.*v)/length(v)) ;end %
function cn=cosndx(arr,ndx); cn=cosangle(arr(ndx+1,:)'-arr(ndx,:)',arr(ndx,:)'-arr(ndx-1,:)') ;end
function v=vectorproduct(x,y);v=[x(2)*y(3)-x(3)*y(2),x(3)*y(1)-x(1)*y(3),x(1)*y(2)-x(2)*y(1)]';end% for 3x1 vectors x and y only
function ndcs=turningpoints(arr);ndcs=zeros(1,5);tp=0; for  i=2:(length(arr)-1);cn=cosndx(arr,i); if cn<0.5; tp=tp+1;ndcs(tp)=i; endif; i=i+1; endfor ;end
function plot_lg(T,y,titlu);figure;subplot(2,2,1);semilogy(T,y);title(titlu);xlabel('time');ylabel('y_1 , y_2 and y_3');subplot(2,2,2);loglog(y(:,1),y(:,2));xlabel('y1'); ylabel('y2');subplot(2,2,3);loglog(y(:,2),y(:,3));xlabel('y2'); ylabel('y3');subplot(2,2,4);loglog(y(:,1),y(:,3));xlabel('y1'); ylabel('y3');end
function plot_it(T,y,titlu);figure;subplot(2,2,1);plot(T,y);title(titlu);xlabel('time');ylabel('y_1 , y_2 and y_3');subplot(2,2,2);plot(y(:,1),y(:,2));xlabel('y1'); ylabel('y2');subplot(2,2,3);plot(y(:,2),y(:,3));xlabel('y2'); ylabel('y3');subplot(2,2,4);plot(y(:,1),y(:,3));xlabel('y1'); ylabel('y3');end
function x=normalize(y);n=1.0/sqrt(sum(y.*y)); x=y*n; end
