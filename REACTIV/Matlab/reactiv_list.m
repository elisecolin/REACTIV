%  function C=reactiv_list(Rep,List,threshold,L,datenum)
%
%  Input:
%  Rep          ::  Repertory where data are.
%  List         ::  List of files containing images in the time serie
%  threshold    ::  Threshold for amplitude - Optional : let it to 0 if you
%  do not want to specify
%  L            ::  equivalent number of looks - Optional : let it to 0 if you do not want to specify 
%               ::  If not specified, estimated by default from time-series
%  datenum      ::  Vector of dates in datenum matlab format. Optional. If
%  not specified, dates are equally distributed over time and images are supposed to
%  be ordered.
%
%  
%  (c) Elise Koeniguer 
%  


function C=reactiv_list(Rep,List,threshold,L,dateN);
fprintf('--------------------------------------------- \n')
fprintf('REACTIV ------------------------------------- \n')

    N=length(List);
    X=imread([Rep,List(1).name]);
    [nx,ny]=size(X);
   
    if (nargin<3 || threshold==0)
        threshold=5*mean(X(:)); 
        fprintf('Automatic Threshold over Amplitude: %f \n',threshold)
    else
        fprintf('Fixed Threshold over Amplitude: %f \n',threshold)
    end
    
    

    
    
    if nargin<5
        dateN=(1:N); % If no dates are defined, there are assumed equally distributed
        str1='Date Begin';
        str2='Date End';
        fprintf('Dates are assumed to be equally distributed \n')
    else
        str1=datestr(min(dateN),26);
        str2=datestr(max(dateN),26);
        fprintf('Dates defined between: %s and %s \n',str1,str2)
        
    end
    ds=max(dateN)-min(dateN);
    
    
    m1=X;
    m2=X.^2;
    Imax=X;
    Kmax=zeros(nx,ny);

    
    for k=2:N
        X=imread([Rep,List(k).name]);
        m1=m1+X;
        m2=m2+X.^2;
        Matrix_indicek=(dateN(k)-min(dateN))*ones(nx,ny);
        Kmax(X>Imax)=Matrix_indicek(X>Imax); % Replace Indice where X is higher than current value
        Imax=max(Imax,X);
    end
    m1=m1/N;
    m2=m2/N;

    R=sqrt(m2-m1.^2)./m1;
    gam=mean(R(:));
    
    
    if (nargin<4||L==0)
        % Simplified estimation of L trough coefficient of variation for Rayleigh-Nakagami
        L= ((0.991936+0.067646*gam-0.098888*gam^2 -0.048320  *gam^3)/(0.001224-0.034323*gam+4.305577*gam^2-1.163498*gam^3));
        fprintf('Automatic Estimation of L, Equivalent Number of Looks: %f \n',L)
    end
    
    fprintf('Normalization of Coefficient of Variation and Intensity \n')
    % R Normalization 
    CV=sqrt((L.*gamma(L).^2./(gamma(L+0.5).^2))-1); % theretical mean value
    num=(L.*gamma(L).^4.*(4*L.^2.*gamma(L).^2-4*L.*gamma(L+1/2).^2-gamma(L+1/2).^2));
    den=(gamma(L+1/2).^4.*(L.*gamma(L).^2-gamma(L+1/2).^2));
    alpha=1/4*num./den; % theretical standard deviation value
    
   
    R=(R-CV)/(alpha/sqrt(N))/10+0.25;
    R=(R>1)+(R<1).*R;   % Cast Coefficient of Varation R max to 1.
    R=(R>0).*R;
      
    I=(Imax/threshold); % Cast Intensity to threshold. 
    I=(I<1).*I+(I>1);

    hue=Kmax/ds;


    HSV(:,:,1)=hue;
    HSV(:,:,2)=R;
    HSV(:,:,3)=I;
    fprintf('Conversion from HSV to RGB \n')

    C=hsv2rgb(HSV);

    
    figure;
    imagesc(C)
    axis off
    h=colorbar('XTickLabel',{str1,str2},'location','SouthOutside');
    colormap(hsv)
    set(h,'XTick',[0,1]);
    set(h,'FontSize',14)
    

