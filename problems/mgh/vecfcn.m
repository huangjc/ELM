function fvec = vecfcn( x, nprob, m )
% f: R^n --> R^m
n = length( x );

switch nprob
    % case 1-35 are examples from more81 
        
    case 1 % Rosenbrock function (more81-1)
        fvec = zeros(2,1);
        fvec(1) = 10*(x(2)-x(1)^2);
        fvec(2) = 1-x(1);        
    
    case 2 % Freudenstein and Roth function 
        fvec = zeros(2,1);
        fvec(1) = -13 + x(1) + ((5-x(2))*x(2)-2)*x(2);
        fvec(2) = -29 + x(1) + ((x(2)+1)*x(2)-14)*x(2);
    
    case 3 % Powell badly scaled function (more81-3)
        fvec = zeros(2,1);
        fvec(1) = (10^4)*x(1)*x(2)-1;
        fvec(2) = exp(-x(1))+exp(-x(2))-1.0001;
        
    case 4 % Brown badly scaled function 
        fvec = zeros(3,1);
        fvec(1) = x(1)-10^6;
        fvec(2) = x(2)-2*10^(-6);
        fvec(3) = x(1)*x(2)-2;
        
    case 5 % Beale function
        fvec = zeros(3,1);
        fvec(1) = 1.5 - x(1)*(1-x(2));
        fvec(2) = 2.25 - x(1)*(1-x(2)^2);
        fvec(3) = 2.625 - x(1)*(1-x(2)^3);
        
    case 6 % Jenrich and Sampson function
        fvec = zeros(m,1);
        for i = 1:m;
            fvec(i) = 2 + 2*i -(exp(i*x(1)) + exp(i*x(2)));
        end
        
    case 7 %  Helical valley function (more81-7)
        fvec = zeros(3,1);
        if (x(1)>0)
            temp1=atan(x(2)/x(1))/2/pi;
        elseif (x(1)<0) 
            temp1=atan(x(2)/x(1))/2/pi+0.5;
        else % x(1) == 0 not defined?
            temp1 = 1/4;
        end
        temp2=sqrt(x(1)^2+x(2)^2);
        fvec(1)=10*(x(3)-10*temp1);
        fvec(2)=10*(temp2-1);
        fvec(3)=x(3);
        
    case 8 % Bard function
        fvec = zeros(15,1);
        y = [0.14; 0.18; 0.22; 0.25; 0.29; 0.32; 0.35; 0.39;...
            0.37; 0.58; 0.73; 0.96; 1.34; 2.10; 4.39 ];
        u = [1:15];
        v = 16 - [1:15];
        w = min([u;v]);
        for i = 1:15;
            fvec(i) = y(i) - ( x(1) + u(i)/(v(i)*x(2)+w(i)*x(3)) );
        end
        
    case 9 % Gaussian function
        fvec = zeros(15,1);
        y = [0.0009; 0.0044; 0.0175; 0.0540; 0.1295; 0.2420; 0.3521; 0.3989];
        y = [ y; y(7:-1:1)];
        t = ( 8 - [1:15] )/2;
        for i = 1:15;
            fvec(i) = x(1)*exp(-x(2)*(t(i)-x(3))^2 /2 ) -y(i);
        end
        
    case 10 % Meyer function
        fvec = zeros(16,1);
        t = 45 + 5*[1:16];
        y = [34780; 28610; 23650; 19630; 16370; 13720; 11540; 9744;
            8261; 7030; 6005; 5147; 4427; 3820; 3307; 2872 ];
        for i = 1:16; 
            fvec(i) = x(1)*exp(x(2)/(t(i)+x(3))) -y(i);
        end
        
    case 11 % Gulf research and development function
        fvec = zeros(m,1);
        t = [1:m]/100;
        y = 25+(-50*log(t)).^(2/3);
        for i = 1:m;
            fvec(i) = exp( -abs( y(i)*m*i*x(2) )^(x(3))/x(1) ) - t(i);
        end
        
    case 12 % Box three-dimensional function
        fvec = zeros(m,1);
        t = 0.1*[1:m];
        for i = 1:m;
            fvec(i) = exp(-t(i)*x(1)) - exp(-t(i)*x(2)) - x(3)*(exp(-t(i))-exp(-10*t(i)));
        end
        
    case 13 %  Powell singular function (more81-13)
        fvec = zeros(4,1);
        fvec(1)=x(1)+10*x(2);
        fvec(2)=sqrt(5)*(x(3)-x(4));
        fvec(3)=(x(2)-2*x(3))^2;
        fvec(4)=sqrt(10)*(x(1)-x(4))^2;
        
    case 14 %  Wood function (more81-14)
        fvec = zeros(6,1);
        temp1=x(2)-x(1)^2;
        temp2=x(4)-x(3)^2;
        fvec(1)=10*temp1;
        fvec(2)=1-x(1);
        fvec(3)=sqrt(90)*temp2;
        fvec(4)=1-x(3);
        fvec(5)=sqrt(10)*(x(2)+x(4)-2);
        fvec(6)=(x(2)-x(4))/sqrt(10);

    case 15 % Kowalik and Osborne function
        y = [0.1957; 0.1947; 0.1735; 0.1600; 0.0844; 0.0627;
            0.0456; 0.0342; 0.0323; 0.0235; 0.0246];
        u = [4; 2; 1; 0.5; 0.25; 0.167; 0.125; 0.1; 0.0833; 0.0714; 0.0625];
        fvec = zeros(11,1);
        for i = 1:11;
            fvec(i) = y(i) - x(1)*(u(i)^2 +u(i)*x(2)) / (u(i)^2+u(i)*x(3)+x(4));
        end
      
    case 16 % Brown and Dennis function
        t = [1:m]/5;
        fvec = zeros(m,1);
        for i = 1:m;
            fvec(i) = (x(1)+t(i)*x(2)-exp(t(i)))^2 + (x(3)+x(4)*sin(t(i))-cos(t(i)))^2;
        end
        
    case 17 % Osborne 1 function
        t = 10*([1:33]-1);
        y = [0.844; 0.908; 0.932; 0.936; 0.925; 0.908; 0.881; 0.850; 0.818; 
            0.784; 0.751; 0.718; 0.685; 0.658; 0.628; 0.603; 0.580; 0.558;
            0.538; 0.522; 0.506; 0.490; 0.478; 0.467; 0.457; 0.448; 0.438; 
            0.431; 0.424; 0.420; 0.414; 0.411; 0.406];
        fvec = zeros(33,1); 
        for i = 1:33;
            fvec(i) = y(i) - (x(1)+x(2)*exp(-t(i)*x(4))+x(3)*exp(-t(i)*x(5)));
        end
        
    case 18 % Biggs EXP6 function
        t = 0.1*[1:m];
        fvec = zeros(m,1);
        for i = 1:m;  
            yi = exp(-t(i)) - 5*exp(-10*t(i)) + 3*exp(-4*t(i));
            fvec(i) = x(3)*exp(-t(i)*x(1)) - x(4)*exp(-t(i)*x(2)) +x(6)*exp(-t(i)*x(5)) -yi;
        end
        
    case 19 % Osborne 2 function 
        t = ([1:m]-1)/10;
        y = [1.366; 1.191; 1.112; 1.013; 0.991; 0.885; 0.831; 0.847; 0.786; 0.725; 0.746; 0.679; 0.608; 0.655; 0.616; 0.606; 0.602; 0.626; 0.651; 0.724; 0.649; 0.649;
            0.694; 0.644; 0.624; 0.661; 0.612; 0.558; 0.533; 0.496; 0.500; 0.423; 0.395; 0.375; 0.372; 0.391; 0.396; 0.405; 0.428; 0.429; 0.523; 0.562; 0.607; 0.653; 
            0.672; 0.708; 0.633; 0.668; 0.645; 0.632; 0.591; 0.559; 0.597; 0.625; 0.739; 0.710; 0.729; 0.720; 0.636; 0.581; 0.428; 0.292; 0.162; 0.098; 0.054 ];
        fvec  = zeros(m,1);
        for i = 1:m; 
            fvec(i) = y(i) - (x(1)*exp(-t(i)*x(5)) + x(2)*exp(-(t(i)-x(9))^2*x(6)) + x(3)*exp(-(t(i)-x(10))^2*x(7)) + x(4)*exp(-(t(i)-x(11))^2*x(8)));
        end
        
    case 20 %  Watson function (more81-20)
        fvec = zeros(m,1);
        t = [1:29]'/29;
        for i=1:29
            t1 = t(i).^[0:n-1]';
            fvec(i) = sum( [1:n-1]'.*x(2:n).*t1(1:n-1) ) - sum( x.*t1 )^2 - 1;
%             ti=i/29;
%             sum1=0;
%             temp=1;
%             for j=2:n
%                 sum1=sum1+(j-1)*temp*x(j);
%                 temp=ti*temp;
%             end
%             sum2=0;
%             temp=1;
%             for j=1:n
%                 sum2=sum2+temp*x(j);
%                 temp=ti*temp;
%             end
%             fvec1=sum1-sum2^2-1;
%             if fvec(i) - fvec1 == 0
%                 error('')
%             end
        end
        fvec(30) = x(1);
        fvec(31) = x(2)-x(1)^2-1;
        
    case 21 % Extended Rosenbrock function
        % Not sure what x_(i-1) is when i = 1, assume it is equivalent to
        % Rosenbrock function when n=2
        % so f(2*i-1) = 10*( x(2*i) - x(2*i-1)^2 )?
        fvec = zeros(m,1);
        for i = 1:n
            if mod(i,2) == 1
                fvec(i) = 10*( x(i+1) - x(i)^2 );
            else
                fvec(i) = 1 - x(i-1);
            end
        end
        
    case 22 % Extended Powell singular function
        fvec = zeros(m,1);
        for i = 1:(m/4);
            fvec(4*i-3) = x(4*i-3) + 10*x(4*i-2);
            fvec(4*i-2) = sqrt(5)*(x(4*i-1)-x(4*i));
            fvec(4*i-1) = (x(4*i-2)-2*x(4*i-1))^2;
            fvec(4*i) = sqrt(10)*(x(4*i-3)-x(4*i))^2;
        end
        
    case 23 % Penalty function I
        % m = n+1;
        a = 10^(-5);
        sqrta  = sqrt(a);
        fvec = zeros(m,1);
        for i = 1:n;
            fvec(i) = sqrta*(x(i)-1);
        end
        fvec(n+1) = sum(x.^2)-1/4;
        
    case 24 % Penalty function II
        % m = 2*n;
        a = 10^(-5);
        sqrta  = sqrt(a);
        fvec = zeros(m,1);
        fvec(1) = x(1) - 0.2; 
        for i = 2:n;
            yi = exp(i/10) + exp((i-1)/10);
            fvec(i) = sqrta*( exp(x(i)/10) + exp(x(i-1)/10) - yi );
        end
        for i = (n+1):(2*n-1)
            fvec(i) = sqrta*( exp(x(i-n+1)/10) + exp(-1/10) );
        end
        fvec(2*n) = sum( (n-[1:n]'+1).*(x.^2) ) -1;
         
    case 25 %  Variably dimensioned function (more81-25)
        % m = n+2;
        fvec = zeros(m,1);
        for i=1:n
            fvec(i)=x(i)-1;
        end
        sum1=0;
        for j=1:n
            sum1=sum1+j*(x(j)-1);
        end
        fvec(m-1)=sum1;
        fvec(m)=sum1^2;
        
    case 26 %  Trigonometric function (more81-26)
        sum1=0;
        fvec = zeros(n,1);
        for j=1:n
            fvec(j)=cos(x(j));
            sum1=sum1+fvec(j);
        end
        for k=1:n
            fvec(k)=(n+k)-sin(x(k))-sum1-k*fvec(k);
        end
        
    case 27 %  Brown almost-linear function (more81-27)
        fvec = zeros(n,1);
        sum1=-(n+1);
        prod=1;
        for j=1:n
            sum1=sum1+x(j);
            prod=x(j)*prod;
        end
        for k=1:n
            fvec(k)=x(k)+sum1;
        end
        fvec(n)=prod-1;

    case 28 %  Discrete boundary value function (more81-28)
        fvec = zeros(n,1);
        h=1/(n+1);
        for k=1:n
            temp=(x(k)+k*h+1)^3;
            temp1=0;
            if k ~=1
                temp1=x(k-1);
            end
            temp2=0;
            if k ~=n
                temp2=x(k+1);
            end
            fvec(k)=2*x(k)-temp1-temp2+temp*(h^2)/2;
        end

    case 29 %  Discrete integral equation function (more81-29)
        fvec = zeros(n,1);
        h=1/(n+1);
        for k=1:n
            tk=k*h;
            sum1=0;
            for j=1:k
                tj=j*h;
                temp=(x(j)+tj+1)^3;
                sum1=sum1+tj*temp;
            end
            sum2=0;
            kp1=k+1;
            if n >= kp1
                for j=kp1:n
                    tj=j*h;
                    temp=(x(j)+tj+1)^3;
                    sum2=sum2+(1-tj)*temp;
                end
            end
            fvec(k)=x(k)+h*((1-tk)*sum1+tk*sum2)/2;
        end

    %  Broyden tridiagonal function (more81-30)
    case 30
        fvec = zeros(n,1);
        for k=1:n
            temp=(3-2*x(k))*x(k);
            temp1=0;
            if k ~= 1
                temp1=x(k-1);
            end
            temp2=0;
            if k ~= n
                temp2=x(k+1);
            end
            fvec(k)=temp-temp1-2*temp2+1;
        end

    %  Broyden banded function (more81-31)
    case 31
        fvec = zeros(n,1);
        ml=5;
        mu=1;
        for k=1:n
            k1=max([1;k-ml]);
            k2=min([k+mu;n]);
            temp=0;
            for j=k1:k2
                if j ~=k
                    temp=temp+x(j)*(1+x(j));
                end
            end
            fvec(k)=x(k)*(2+5*x(k)^2)+1-temp;
        end
        
    case 32 % Linear function - full rank
        fvec = zeros(m,1);
        sumx = sum(x);
        for i = 1:n;
            fvec(i) = x(i) - 2/m*sumx -1;
        end
        for i = (n+1):m;
            fvec(i) = -2/m*sumx -1;
        end
        
    case 33 % Linear function - rank 1
        fvec = zeros(m,1);
        for i = 1:m;
            fvec(i) = i*sum([1:n]'.*x) -1;
        end
        
    case 34 % Linear function - rank 1 with 0 columns and rows
        fvec = zeros(m,1);
        fvec(1) = -1;
        fvec(m) = -1;
        for i = 2:(m-1);
            fvec(i) = (i-1)*sum([2:(n-1)]'.*x(2:(n-1))) -1;
        end
        
    case 35 % Chebyquad function (more81-35)
        fvec = zeros(n,1);
        for j=1:n
            temp1=1;
            temp2=2*x(j)-1;
            temp=2*temp2;
            for i=1:n
                fvec(i)=fvec(i)+temp2;
                ti=temp*temp2-temp1;
                temp1=temp2;
                temp2=ti;
            end
        end
        tk=1/n;
        iev=-1;
        for k=1:n
            fvec(k)=tk*fvec(k);
            if iev>0
                fvec(k)=fvec(k)+1/(k^2-1);
            end
            iev=-iev;
        end
        
%     case 15
%         temp=x(1)^2+x(2)^2;
%         f(1)=-1-1/(temp+1);
%         f(2)=sin(x(1))+cos(x(2))-3;
% 
%     case 16
%         temp=x(1)^2+x(2)^2;
%         f(1)=1-1/(temp+1);
%         f(2)=sin(x(1))+cos(x(2))-1;
% 
%     case 17
%         %     if x(1)^2+(x(2)+1)^2<=1/4
%         if x(1)^2+(x(2)+1)^2<=1
%             temp=x(1)^2+x(2)^2;
%             f(1)=1-1/(temp+1);
%             f(2)=sin(x(1))+cos(x(2))-1;
%         else
%             temp=x(1)^2+x(2)^2;
%             f(1)=-1-1/(temp+1);
%             f(2)=sin(x(1))+cos(x(2))-3;
%         end
% 
%     case 18
%         temp=x(1)^2+20*x(2)^2;
%         f(1)=-1-1/(temp+1);
%         f(2)=sin(x(1))+cos(x(2))-3;
% 
%     case 19
%         f(1)=x(1)^2+x(2)^2+1;
% 
%     case 20
%         temp=x(1)^2+x(2)^2;
%         f(1)=-1-1/(temp+1);
%         f(2)=sin(x(1))+cos(x(2))-3;
% 
%         temp2=(x(1)-1.2)^2+x(2)^2-0.4^2;
%         if temp2<0
%             gvec(1)=-350*x(1)+355;
%             gvec(2)=50*x(2)-45;
%             f(1)=f(1)-temp2^2*gvec(1);
%             f(2)=f(2)-temp2^2*gvec(2);
%         end
% 
%      case 21
%         temp=x(1)^2+x(2)^2;
%         f(1)=-1-1/(temp+1);
%         f(2)=sin(x(1))+cos(x(2))-3;
% 
%         temp2=(x(1)-1.2)^2+x(2)^2-0.4^2;
%         k=1/0.4^2+10;
%         if temp2<0
%             f(1)=f(1)*(1+k*temp2);
%             f(2)=f(2)*(1+k*temp2);
%         end   
% 
%      case 22    
%         f(1)=x(2)^2;
%         f(2)=x(2)*(1+x(1)^2);
% 
%      case 23    
%         f(1)=x(1)^2+x(2);
%         f(2)=x(2)^2;    
end

% if size( fvec, 2 ) ~= 1;
%     fvec = fvec';
% end