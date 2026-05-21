function setup = SetupSigns(nvar,m,nlag,horizons,h)


NS        = 1 + numel(horizons);

% setup
setup.nvar=nvar;
setup.m=m;
setup.h=h;
setup.e = eye(nvar);
setup.nlag     = nlag;
setup.horizons = horizons;
setup.nex      = m-nvar*nlag;

%==========================================================================
%% identification: 
% L_0 is the matrix of contemporaneous impulse responses
% L_0(i,j) is the response of variable i to shock j
% L_0(i,j) S_{i,j} L_0 e_j 
%==========================================================================

setup.S = cell([setup.nvar,1]);
for ii=1:nvar
    setup.S{ii}=zeros(0,setup.nvar*NS);
end

% shock 1 = demand shock
nr1=14; % number of restrictions imposed on shock 1   
setup.Ss{1} = zeros(nr1,setup.nvar);
setup.Ss{1}(1,1)=1;   % GDP + 
setup.Ss{1}(2,10)=1;  % GDP DEF +
setup.Ss{1}(3,11)=1;  % PCE INDEX + 
setup.Ss{1}(4,12)=1;  % PCEXFE INDEX +
setup.Ss{1}(5,13)=1;  % CPI +
setup.Ss{1}(6,14)=1;  % CPIXFE +
setup.Ss{1}(7,19)=-1; % UNEMPLOYMENT -
setup.Ss{1}(8,20)=1;  % INDPRO +
setup.Ss{1}(9,21)=1;  % CAP UTIL +
setup.Ss{1}(10,25)=1; % FED FUNDS +
setup.Ss{1}(11,26)=1; % 3-month t-bill +
setup.Ss{1}(12,30)=1; % PRIME RATE +
setup.Ss{1}(13,4)=-1;setup.Ss{1}(13,1)=+1; % NON RESIDENTIAL minus GDP
setup.Ss{1}(14,7)=-1;setup.Ss{1}(14,1)=+1; % GOV SPENDING minus GDP


% shock 2 = investment shock
nr2=15;
setup.Ss{2}= zeros(nr2,nvar);
setup.Ss{2}(1,1)=1;    % GDP + 
setup.Ss{2}(2,10)=1;   % GDP DEF +
setup.Ss{2}(3,11)=1;   % PCE index +
setup.Ss{2}(4,12)=1;   % PCEXFE index +
setup.Ss{2}(5,13)=1;   % CPI +
setup.Ss{2}(6,14)=1;   % CPIXFE +
setup.Ss{2}(7,19)=-1;  % U -
setup.Ss{2}(8,20)=1;   % INDPRO +
setup.Ss{2}(9,21)=1;   % CAP UTIL +
setup.Ss{2}(10,25)=1;  % FED FUNDS + 
setup.Ss{2}(11,26)=1;  % 3-MONTH T-BILL +
setup.Ss{2}(12,30)=1;  % PRIME RATE + 
setup.Ss{2}(13,34)=-1;  % S&P 500 - 
setup.Ss{2}(14,4)=1;setup.Ss{2}(14,1)=-1;
setup.Ss{2}(15,7)=-1;setup.Ss{2}(15,1)=+1;

% shock 3 = financial shock
nr3=15;
setup.Ss{3}= zeros(nr3,nvar);
setup.Ss{3}(1,1)=1;    % GDP + 
setup.Ss{3}(2,10)=1;   % GDP DEF +
setup.Ss{3}(3,11)=1;   % PCE index +
setup.Ss{3}(4,12)=1;   % PCEXFE index +
setup.Ss{3}(5,13)=1;   % CPI +
setup.Ss{3}(6,14)=1;   % CPIXFE +
setup.Ss{3}(7,19)=-1;  % U -
setup.Ss{3}(8,20)=1;   % INDPRO +
setup.Ss{3}(9,21)=1;   % CAP UTIL +
setup.Ss{3}(10,25)=1;  % FED FUNDS + 
setup.Ss{3}(11,26)=1;  % 3-MONTH T-BILL +
setup.Ss{3}(12,30)=1;  % PRIME RATE + 
setup.Ss{3}(13,34)=1;  % S&P 500 +
setup.Ss{3}(14,4)=1; setup.Ss{3}(14,1)=-1;
setup.Ss{3}(15,7)=-1;setup.Ss{3}(15,1)=+1;


% shock 4 = monetary policy shock
nr4=19;
setup.Ss{4}= zeros(nr4,nvar);
setup.Ss{4}(1,1)=-1;    % GDP -
setup.Ss{4}(2,10)=-1;   % GDP DEF -
setup.Ss{4}(3,11)=-1;   % PCE index -
setup.Ss{4}(4,12)=-1;   % PCEXFE index -
setup.Ss{4}(5,13)=-1;   % CPI -
setup.Ss{4}(6,14)=-1;   % CPIXFE -
setup.Ss{4}(7,18)=-1;   % EMPL -
setup.Ss{4}(8,19)= 1;   % U +
setup.Ss{4}(9,20)=-1;   % INDPRO -
setup.Ss{4}(10,21)=-1;  % CAP UTIL -
setup.Ss{4}(11,25)=1;   % FED FUNDS + 
setup.Ss{4}(12,26)=1;   % 3-MONTH T-BILL +
setup.Ss{4}(13,27)=1;   % 2-YEAR T-NOTE +
setup.Ss{4}(14,28)=1;   % 5-YEAR T-NOTE +
setup.Ss{4}(15,29)=1;   % 10-YEAR T-NOTE +
setup.Ss{4}(16,30)=1;   % PRIME RATE + 
setup.Ss{4}(17,31)=1;   % AAA CORPORTATe + 
setup.Ss{4}(18,32)=1;   % BAA CORPORTATe + 
setup.Ss{4}(19,34)=-1;  % S&P 500 +

% shock 5 = government spending shock
nr5=14;
setup.Ss{5}= zeros(nr5,nvar);
setup.Ss{5}(1,1)= 1;   % GDP +
setup.Ss{5}(2,7)= 1;   % GOV SPENDING +
setup.Ss{5}(3,8)= -1;  % FED BUDGET +
setup.Ss{5}(4,9)= 1;   % FED TAX RECEIPTS +
setup.Ss{5}(5,10)=1;   % GDP DEF +
setup.Ss{5}(6,11)=1;   % PCE index +
setup.Ss{5}(7,12)=1;   % PCEXFE index +
setup.Ss{5}(8,13)=1;   % CPI +
setup.Ss{5}(9,14)=1;   % CPIXFE +
setup.Ss{5}(10,19)=-1; % UNEMPLOYMENT -
setup.Ss{5}(11,25)=1;  % FFR +
setup.Ss{5}(12,26)=1;  % 3-MONTH T-BILL + 
setup.Ss{5}(13,30)=1;  % PRIME +
setup.Ss{5}(14,7)=1;setup.Ss{5}(14,1)=-1;


% shock 6 = technology shock
nr6=12;
setup.Ss{6}= zeros(nr6,nvar);
setup.Ss{6}(1,1)= 1;  % GDP +
setup.Ss{6}(2,2)= 1;  % CONS EXPENDITURE +
setup.Ss{6}(3,4)= 1;  % NON RESIDENTIAL +
setup.Ss{6}(4,10)=-1; % GDP DEF -
setup.Ss{6}(5,11)=-1;    % PCE index -
setup.Ss{6}(6,12)=-1;    % PCEXFE index -
setup.Ss{6}(7,13)=-1;    % CPI -
setup.Ss{6}(8,14)=-1;    % CPIXFE -
setup.Ss{6}(9,15)= 1;    % HOURLY WAGE +
setup.Ss{6}(10,16)= 1;   % LABOR PROD +
setup.Ss{6}(11,17)= 1;   % UTILIZATION ADJUSTED +
setup.Ss{6}(12,19)=-1;   % UNEMPLOYMENT -


% shock 7 = labor supply shock
nr7=8;
setup.Ss{7}= zeros(nr7,nvar);
setup.Ss{7}(1,1)= 1;    % GDP +
setup.Ss{7}(2,10)=-1;   % GDP DEF -
setup.Ss{7}(3,11)=-1;   % PCE index -
setup.Ss{7}(4,12)=-1;   % PCEXFE index -
setup.Ss{7}(5,13)=-1;   % CPI -
setup.Ss{7}(6,14)=-1;   % CPIXFE -
setup.Ss{7}(7,15)= -1;  % HOURLY WAGE -
setup.Ss{7}(8,19)= 1;   % UNEMPLOYMENT +


% shock 8= wage bargaining shock
nr8=8;
setup.Ss{8}= zeros(nr8,nvar);
setup.Ss{8}(1,1)=  1;   % GDP +
setup.Ss{8}(2,10)=-1;   % GDP DEF -
setup.Ss{8}(3,11)=-1;   % PCE index -
setup.Ss{8}(4,12)=-1;   % PCEXFE index -
setup.Ss{8}(5,13)=-1;   % CPI -
setup.Ss{8}(6,14)=-1;   % CPIXFE -
setup.Ss{8}(7,15)= -1;   % HOURLY WAGE -
setup.Ss{8}(8,19)= -1; % UNEMPLOYMENT -

% shock 9= oil
nr9=13;
setup.Ss{9}= zeros(nr9,nvar);
setup.Ss{9}(1,1)= 1; % GDP +
setup.Ss{9}(2,2)= 1; % CONS EXPENDITURE +
setup.Ss{9}(3,4)= 1; % NON RESIDENTIAL +
setup.Ss{9}(4,10)=-1;   % GDP DEF -
setup.Ss{9}(5,11)=-1;   % PCE index -
setup.Ss{9}(6,12)=-1;   % PCEXFE index -
setup.Ss{9}(7,13)=-1;   % CPI -
setup.Ss{9}(8,14)=-1;   % CPIXFE -
setup.Ss{9}(9,15)= 1;   % HOURLY WAGE +
setup.Ss{9}(10,16)= 1;   % LABOR PROD +
setup.Ss{9}(11,17)= 1;   % UTILIZATION ADJUSTED +
setup.Ss{9}(12,19)= 1; % UNEMPLOYMENT 
setup.Ss{9}(13,35)= -1; % oil price 


% shock 10= sentiment
nr1=11;
setup.Ss{10}= zeros(nr1,nvar);
setup.Ss{10}(1,1)=1; % GDP + 
setup.Ss{10}(2,2)=1; % PCE + 
setup.Ss{10}(3,3)=1; % residential + 
setup.Ss{10}(4,4)=1; % non-residential + 
setup.Ss{10}(5,10)=1; % GDP DEF +
setup.Ss{10}(6,11)=1; % PCE INDEX + 
setup.Ss{10}(7,12)=1; % PCEXFE INDEX +
setup.Ss{10}(8,13)=1; % CPI +
setup.Ss{10}(9,14)=1; % CPIXFE +
setup.Ss{10}(10,19)=1; % UNEMPLOYMENT -
setup.Ss{10}(11,34)=1; % s&p +




