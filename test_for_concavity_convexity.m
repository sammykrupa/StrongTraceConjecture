%This script is used for checking conditions of convexity/concavity for the
%pressure





%Check between X(1,1,2) and X(1,1,4)
%Given four  numbers $a,b,c,d\in\mathbb{R}$ with $c>0$  is it possible to find a strictly increasing, strictly concave function $f$ such that $f(0)=a$ and $f(c)=b$ and the integral of $f$ from 0 to $c$ is equal to $d$? Give also the case for $f$ strictly convex
% the function $f$ in this case is actually $-p$ ($p$ is the pressure)
c=X(1,1,4)-X(1,1,2);
a=X(2,2,2);
b=X(2,2,4);
d=(X(3,1,4)-X(2,1,4)*X(2,1,4)*.5)-(X(3,1,2)-X(2,1,2)*X(2,1,2)*.5);

disp("Check between X(1,1,2) and X(1,1,4)")
checkIncreasingShape(a, b, c, d)


%check between X(1,1,4) and X(1,1,5)
%Given four  numbers $a,b,c,d\in\mathbb{R}$ with $c>0$  is it possible to find a strictly increasing, strictly CONVEX function $f$ such that $f(0)=a$ and $f(c)=b$ and the integral of $f$ from 0 to $c$ is equal to $d$? 
% the function $f$ in this case is actually $-p$ ($p$ is the pressure)
c=X(1,1,5)-X(1,1,4);
a=X(2,2,4);
b=X(2,2,5);
d=(X(3,1,5)-X(2,1,5)*X(2,1,5)*.5)-(X(3,1,4)-X(2,1,4)*X(2,1,4)*.5);


disp("check between X(1,1,4) and X(1,1,5)")
checkIncreasingShape(a, b, c, d)


%check between X(1,1,5) and X(1,1,1)
c=X(1,1,1)-X(1,1,5);
a=X(2,2,5);
b=X(2,2,1);
d=(X(3,1,1)-X(2,1,1)*X(2,1,1)*.5)-(X(3,1,5)-X(2,1,5)*X(2,1,5)*.5);
disp("check between X(1,1,5) and X(1,1,1)")
checkIncreasingShape(a, b, c, d)

%check between X(1,1,1) and X(1,1,3)
c=X(1,1,3)-X(1,1,1);
a=X(2,2,1);
b=X(2,2,3);
d=(X(3,1,3)-X(2,1,3)*X(2,1,3)*.5)-(X(3,1,1)-X(2,1,1)*X(2,1,1)*.5);

disp("check between X(1,1,1) and X(1,1,3)")
checkIncreasingShape(a, b, c, d)

%check between X(1,1,3) and X(1,1,6)
c=X(1,1,6)-X(1,1,3);
a=X(2,2,3);
b=X(2,2,6);
d=(X(3,1,6)-X(2,1,6)*X(2,1,6)*.5)-(X(3,1,3)-X(2,1,3)*X(2,1,3)*.5);

disp("check between X(1,1,3) and X(1,1,6)")
checkIncreasingShape(a, b, c, d)