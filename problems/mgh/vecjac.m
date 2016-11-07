function fjac=vecjac( x, nprob, m )

n = length( x );

switch nprob
    case 1 %  Rosenbrock function (more81-1)
        fjac(1,1)=-20*x(1);
        fjac(1,2)=10;
        fjac(2,1)=-1;
        fjac(2,2)=0;
    
    case 2 % Freudenstein and Roth function 
        fjac(1,1)=1;
        fjac(1,2)=-3*x(2)^2+10*x(2)-2;
        fjac(2,1)=1;
        fjac(2,2)=3*x(2)^2+2*x(2)-14;
        
    case 3 %  Powell badly scaled function (more81-3)
        fjac(1,1)=10^4*x(2);
        fjac(1,2)=10^4*x(1);
        fjac(2,1)=-exp(-x(1));
        fjac(2,2)=-exp(-x(2));
        
    case 4 % Brown badly scaled function 
        fjac(1,1)=1;
        fjac(1,2)=0;
        fjac(2,1)=0;
        fjac(2,2)=1;
        fjac(3,1)=1;
        fjac(3,2)=1;
        
    case 5 % Beale function
        fjac(1,1)=x(2)-1;
        fjac(1,2)=x(1);
        fjac(2,1)=x(2)^2-1;
        fjac(2,2)=2*x(1)*x(2);
        fjac(3,1)=x(2)^3-1;
        fjac(3,2)=3*x(1)*x(2)^2;
        
    case 6 % Jenrich and Sampson function
        fjac = zeros(m,n);
        for i = 1:m
            fjac(i,1)=-i*exp(i*x(1));
            fjac(i,2)=-i*exp(i*x(2));
        end
        
    case 7 %  Helical valley function (more81-7)
        tpi=8*atan(1);
        temp=x(1)^2+x(2)^2;
        temp1=tpi*temp;
        temp2=sqrt(temp);
        fjac(1,1)=100*x(2)/temp1;
        fjac(1,2)=-100*x(1)/temp1 ;
        fjac(1,3)=10;
        fjac(2,1)=10*x(1)/temp2;
        fjac(2,2)=10*x(2)/temp2 ;
        fjac(2,3)=0;
        fjac(3,1)=0;
        fjac(3,2)=0;
        fjac(3,3)=1;
        
    case 8 % Bard function
        fjac = zeros(m,3);
        u = [1:m];
        v = 16 - [1:m];
        w = min([u;v]);
        for i = 1:m;
            fjac(i,1) = -1;
            fjac(i,2) = u(i)*v(i)/(v(i)*x(2)+w(i)*x(3))^2;
            fjac(i,3) = u(i)*w(i)/(v(i)*x(2)+w(i)*x(3))^2;
        end
        
    case 9 % Gaussian function
        fjac = zeros(m,3);
        t = ( 8 - [1:m] )/2;
        for i = 1:m;
            fjac(i,1) = exp(-x(2)*(t(i)-x(3))^2 /2 );
            fjac(i,2) = -x(1)*exp(-x(2)*(t(i)-x(3))^2 /2 ) * (t(i)-x(3))^2 /2 ;
            fjac(i,3) = -x(1)*exp(-x(2)*(t(i)-x(3))^2 /2 ) * x(2)*(t(i)-x(3));
        end
        
    case 10 % Meyer function
        fjac = zeros(m,3);
        t = 45 + 5*[1:m];
        for i = 1:m; 
            fjac(i,1) = exp(x(2)/(t(i)+x(3)));
            fjac(i,2) = x(1)*exp(x(2)/(t(i)+x(3))) /(t(i)+x(3));
            fjac(i,3) = -x(1)*exp(x(2)/(t(i)+x(3))) *x(2)/(t(i)+x(3))^2;
        end
        
    case 11 % Gulf research and development function
        fjac = zeros(m,n);
        t = [1:m]/100;
        y = 25+(-50*log(t)).^(2/3);
        for i = 1:m; 
            fjac(i,1) = exp( -abs( y(i)*m*i*x(2) )^(x(3))/x(1) ) * abs( y(i)*m*i*x(2) )^(x(3))/x(1)^2;
            fjac(i,2) = -exp( -abs( y(i)*m*i*x(2) )^(x(3))/x(1) ) * (x(3)/x(2))*abs(y(i)*m*i*x(2) )^(x(3)) /x(1);
            fjac(i,3) = -exp( -abs( y(i)*m*i*x(2) )^(x(3))/x(1) ) * abs( y(i)*m*i*x(2) )^(x(3))*log(abs( y(i)*m*i*x(2) ))/x(1);
        end
        
    case 12 % Box three-dimensional function
        fjac = zeros(m,n);
        t = 0.1*[1:m];
        for i = 1:m;
            fjac(i,1) = -t(i)*exp(-t(i)*x(1));
            fjac(i,2) = t(i)*exp(-t(i)*x(2));
            fjac(i,3) = -(exp(-t(i))-exp(-10*t(i)));
        end
        
    case 13 %  Powell singular function (more81-13)
        fjac = zeros(m,n);
        fjac(1,1)=1;
        fjac(1,2)=10;
        fjac(2,3)=sqrt(5);
        fjac(2,4)=-fjac(2,3);
        fjac(3,2)=2*(x(2)-2*x(3));
        fjac(3,3)=-2*fjac(3,2);
        fjac(4,1)=2*sqrt(10)*(x(1)-x(4));
        fjac(4,4)=-fjac(4,1);

    case 14 %  Wood function (more81-14)
        fjac(1,1)=-20*x(1);
        fjac(1,2)=10;
        fjac(2,1)=-1;
        fjac(3,3)=-2*sqrt(90)*x(3);
        fjac(3,4)=sqrt(90);
        fjac(4,3)=-1;
        fjac(5,2)=sqrt(10);
        fjac(5,4)=sqrt(10);
        fjac(6,2)=1/sqrt(10) ;
        fjac(6,4)=-1/sqrt(10) ;

    case 15 % Kowalik and Osborne function
        fjac = zeros(m,n);
        u = [4; 2; 1; 0.5; 0.25; 0.167; 0.125; 0.1; 0.0833; 0.0714; 0.0625];
        for i = 1:m;
            fjac(i,1) = -(u(i)^2 +u(i)*x(2)) / (u(i)^2+u(i)*x(3)+x(4));
            fjac(i,2) = -x(1)*u(i) / (u(i)^2+u(i)*x(3)+x(4));
            fjac(i,3) = u(i)*x(1)*(u(i)^2 +u(i)*x(2)) / (u(i)^2+u(i)*x(3)+x(4))^2;
            fjac(i,4) = x(1)*(u(i)^2 +u(i)*x(2)) / (u(i)^2+u(i)*x(3)+x(4))^2;
        end
      
    case 16 % Brown and Dennis function
        t = [1:m]/5;
        fjac = zeros(m,n);
        for i = 1:m;
            fjac(i,1) = 2*(x(1)+t(i)*x(2)-exp(t(i)));
            fjac(i,2) = 2*t(i)*(x(1)+t(i)*x(2)-exp(t(i)));
            fjac(i,3) = 2*(x(3)+x(4)*sin(t(i))-cos(t(i)));
            fjac(i,4) = 2*sin(t(i))*(x(3)+x(4)*sin(t(i))-cos(t(i)));
        end
        
    case 17 % Osborne 1 function
        t = 10*([1:m]-1);
        fjac = zeros(m,n); 
        for i = 1:33;
            fjac(i,1) = -1;
            fjac(i,2) = -exp(-t(i)*x(4));
            fjac(i,3) = -exp(-t(i)*x(5));
            fjac(i,4) = t(i)*x(2)*exp(-t(i)*x(4));
            fjac(i,5) = t(i)*x(3)*exp(-t(i)*x(5));
        end
        
    case 18 % Biggs EXP6 function
        t = 0.1*[1:m];
        fjac = zeros(m,n); 
        for i = 1:m;  
            fjac(i,1) = -t(i)*x(3)*exp(-t(i)*x(1));
            fjac(i,2) = t(i)*x(4)*exp(-t(i)*x(2));
            fjac(i,3) = exp(-t(i)*x(1));
            fjac(i,4) = - exp(-t(i)*x(2));
            fjac(i,5) = -t(i)*x(6)*exp(-t(i)*x(5));
            fjac(i,6) = exp(-t(i)*x(5));
        end
        
    case 19 % Osborne 2 function 
        m = 65; 
        t = ([1:m]-1)/10;
        fjac = zeros(m,n); 
        for i = 1:m; 
            fjac(i,1) = - exp(-t(i)*x(5));
            fjac(i,2) = - exp(-(t(i)-x(9))^2*x(6));
            fjac(i,3) = - exp(-(t(i)-x(10))^2*x(7));
            fjac(i,4) = - exp(-(t(i)-x(11))^2*x(8));
            fjac(i,5) = t(i)*x(1)*exp(-t(i)*x(5));
            fjac(i,6) = (t(i)-x(9))^2 * x(2)*exp(-(t(i)-x(9))^2*x(6));
            fjac(i,7) = (t(i)-x(10))^2 * x(3)*exp(-(t(i)-x(10))^2*x(7));
            fjac(i,8) = (t(i)-x(11))^2 * x(4)*exp(-(t(i)-x(11))^2*x(8));
            fjac(i,9) = 2*(t(i)-x(9))*x(6) * x(2)*exp(-(t(i)-x(9))^2*x(6));
            fjac(i,10) = 2*(t(i)-x(10))*x(7) * x(3)*exp(-(t(i)-x(10))^2*x(7));
            fjac(i,11) = 2*(t(i)-x(11))*x(8) * x(4)*exp(-(t(i)-x(11))^2*x(8));
        end

    case 20 %  Watson function (more81-20)
        fjac = zeros(m,n); 
        t = [1:29]'/29;
        for i=1:29
            temp = t(i).^[0:n-1]';
            fjac(i,:) = [ 0; [1:n-1]'.*temp(1:n-1) ] - 2*temp*sum( x.*temp );
        end 
        fjac(30,1) = 1;
        fjac(31,1) = -2*x(1);
        fjac(31,2) = 1;
        
    case 21 % Extended Rosenbrock function
        % Not sure what x_(i-1) is when i = 1, assume it is equivalent to
        % Rosenbrock function when n=2
        % so f(2*i-1) = 10*( x(2*i) - x(2*i-1)^2 )?
        fjac = zeros(m,n);
        for i = 1:m;
            if mod(i,2) == 1
                fjac(i,i) = -20*x(i);
                fjac(i,i+1) = 10;
            else
                fjac(i,i-1) = -1;
            end
        end
        
    case 22 % Extended Powell singular function
        fjac = zeros(m,n); 
        for i = 1:(m/4)
            fjac(4*i-3,4*i-3) = 1;
            fjac(4*i-3,4*i-2) = 10;
            fjac(4*i-2,4*i-1) = sqrt(5);
            fjac(4*i-2,4*i) = -sqrt(5);
            fjac(4*i-1,4*i-2) = 2*(x(4*i-2)-2*x(4*i-1));
            fjac(4*i-1,4*i-1) = -4*(x(4*i-2)-2*x(4*i-1));
            fjac(4*i,4*i-3) = 2*sqrt(10)*(x(4*i-3)-x(4*i));
            fjac(4*i,4*i) = -2*sqrt(10)*(x(4*i-3)-x(4*i));
        end
        
    case 23 % Penalty function I
        fjac = zeros(m,n); 
        for i = 1:n
            fjac(i,i) = sqrt(1e-5);
        end
        for j = 1:n
            fjac(m,j) = 2*x(j);
        end        
        
    case 24 % Penalty function II
        fjac = zeros(m,n); 
        sqrta  = sqrt(1e-5);
        fjac(1,1) = 1; 
        for i = 2:n;
            fjac(i,i-1) = sqrta*exp(x(i-1)/10) / 10;
            fjac(i,i) = sqrta*exp(x(i)/10) / 10;
        end
        for i = (n+1):(2*n-1)
            fjac(i,i-n+1) = sqrta*exp(x(i-n+1)/10) / 10;
        end
        for j = 1:n
            fjac(2*n,j) = 2*(n-j+1)*x(j);
        end
         
    case 25 %  Variably dimensi1d function (more81-25)
        fjac = zeros(m,n);
        for i=1:n
            fjac(i,i)=1;
        end
        for j=1:n
            fjac(n+1,j)=j;
        end
        sum1=sum([1:n]'.*(x-1));
        for j=1:n
            fjac(n+2,j)=2*sum1*j;
        end
        
    case 26 %  Trigonometric function (more81-26)
        fjac = zeros(m,n); 
        for j=1:n
            temp=sin(x(j));
            for k=1:n
                fjac(k,j)=temp;
            end
            fjac(j,j)=(j+1)*temp-cos(x(j));
        end
        
    case 27 %  Brown almost-linear function (more81-27)
        product = prod(x);
        fjac = ones(n) + eye(n);
        for j=1:n
            temp=x(j);
            if temp==0
                temp=1;
                product=1;
                for k=1:n
                    if k~=j
                        product=x(k)*product;
                    end
                end
            end
            fjac(n,j)=product/temp;
        end

    case 28 %  Discrete boundary value function (more81-28)
        fjac = zeros(m,n); 
        h=1/(n+1);
        for k=1:n
            temp=3*(x(k)+k*h+1)^2;
            for j=1:n
                fjac(k,j)=0;
            end
            fjac(k,k)=2+temp*h^2/2;
            if k ~=1
                fjac(k,k-1)=-1;
            end
            if k ~=n
                fjac(k,k+1)=-1;
            end
        end

    case 29 %  Discrete integral equation function (more81-29)
        fjac = zeros(m,n); 
        h=1/(n+1);
        for k=1:n
            tk=k*h;
            for j=1:n
                tj=j*h;
                temp=3*(x(j)+tj+1)^2;
                fjac(k,j)=h*min([tj*(1-tk);tk*(1-tj)])*temp/2;
            end
            fjac(k,k)=fjac(k,k)+1;
        end

    case 30 %  Broyden tridiagonal function (more81-30)
        fjac = zeros(m,n); 
        for k=1:n
            fjac(k,k)=3-4*x(k);
            if k ~= 1
                fjac(k,k-1)=-1;
            end
            if k ~= n
                fjac(k,k+1)=-2;
            end
        end

    case 31 %  Broyden banded function (more81-31)
        ml=5;
        mu=1;
        fjac = zeros(m,n); 
        for k=1:n
            k1=max([1;k-ml]);
            k2=min([k+mu;n]);
            for j=k1:k2
                if j ~=k
                    fjac(k,j)=-(1+2*x(j));
                end
            end
            fjac(k,k)=2+15*x(k)^2;
        end
        
    case 32 % Linear function - full rank
        fjac = zeros(m,n); 
        for i = 1:n
            fjac(i,:) = -2/m;
            fjac(i,i) = 1-2/m;
        end
        for i = n+1:m
            fjac(i,:) = -2/m;
        end
        
    case 33 % Linear function - rank 1
        fjac = zeros(m,n); 
        for i = 1:m
            for j = 1:n
                fjac(i,j) = i*j;
            end
        end
        
    case 34 % Linear function - rank 1 with 0 columns and rows
        fjac = zeros(m,n); 
        for i = 2:(m-1)
            for j = 2:(n-1)
                fjac(i,j) = (i-1)*j;
            end
        end
        
    %  Chebyquad function (more81-35)
    case 35 %  Chebyquad function (more81-35)
        fjac = zeros(m,n); 
        tk=1/n;
        for j=1:n
            temp1=1;
            temp2=2*x(j)-1;
            temp=2*temp2;
            temp3=0;
            temp4=2;
            for k=1:n
                fjac(k,j)=tk*temp4;
                ti=4*temp2+temp*temp4-temp3;
                temp3=temp4;
                temp4=ti;
                ti=temp*temp2-temp1;
                temp1=temp2;
                temp2=ti;
            end
        end
        
%     case {15, 16, 17}
%         temp=x(1)^2+x(2)^2+1;
%         fjac(1,1)=2*x(1)/temp^2;
%         fjac(1,2)=2*x(2)/temp^2;
%         fjac(2,1)=cos(x(1));
%         fjac(2,2)=-sin(x(2));
% 
%     case 18
%         temp=x(1)^2+20*x(2)^2+1;
%         fjac(1,1)=2*x(1)/temp^2;
%         fjac(1,2)=40*x(2)/temp^2;
%         fjac(2,1)=cos(x(1));
%         fjac(2,2)=-sin(x(2));
% 
%     case 19
%         fjac(1,1)=2*x(1);
%         fjac(1,2)=2*x(2);
% 
% 
%     case 20
%         temp=x(1)^2+x(2)^2+1;
%         fjac(1,1)=2*x(1)/temp^2;
%         fjac(1,2)=2*x(2)/temp^2;
%         fjac(2,1)=cos(x(1));
%         fjac(2,2)=-sin(x(2));
% 
%         temp2=(x(1)-1.2)^2+x(2)^2-0.4^2;
%         if temp2<0
%             gjac(1,1)=4*temp2*(x(1)-1.2)*(-350*x(1)+355)-350*temp2^2;
%             gjac(1,2)=4*temp2*x(2)*(-350*x(1)+355);
%             gjac(2,1)=4*temp2*(x(1)-1.2)*(50*x(2)-45);
%             gjac(2,2)=4*temp2*x(2)*(50*x(2)-45)+50*temp2^2;
%             fjac=fjac-gjac;
%         end
% 
%     case 21
%         temp0=x(1)^2+x(2)^2;
%         fvec(1)=-1-1/(temp0+1);
%         fvec(2)=sin(x(1))+cos(x(2))-3;
%         temp1=x(1)^2+x(2)^2+1;
%         fjac(1,1)=2*x(1)/temp1^2;
%         fjac(1,2)=2*x(2)/temp1^2;
%         fjac(2,1)=cos(x(1));
%         fjac(2,2)=-sin(x(2));
% 
%         temp2=(x(1)-1.2)^2+x(2)^2-0.4^2;
%         k=1/0.4^2+10;
%         gcn=1+k*temp2;
%         if temp2<0
%             g(1)=k*2*(x(1)-1.2);
%             g(2)=k*2*x(2);
%             fjac=gcn*fjac+[fvec(1);fvec(2)]*g;
%         end
% 
%     case 22
%         fjac(1,1)=0;
%         fjac(1,2)=2*x(2);
%         fjac(2,1)=2*x(1)*x(2);
%         fjac(2,2)=1+x(1)^2;
% 
%     case 23
%         fjac(1,1)=2*x(1);
%         fjac(1,2)=1;
%         fjac(2,1)=0;
%         fjac(2,2)=2*x(2);    

end


