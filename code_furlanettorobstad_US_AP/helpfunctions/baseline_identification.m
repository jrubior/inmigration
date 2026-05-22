% shock 1 = business cycle
nr1=3; % number of restrictions imposed on shock 1   
setup.Ss{1} = zeros(nr1,setup.nvar);
setup.Ss{1}(1,1)=1;   % GDP + 
setup.Ss{1}(2,2)=1;  % Real Wage +
setup.Ss{1}(3,3)=1;  % Labor Force Participation + 



% shock 2 = wage bargaining
nr2=3;
setup.Ss{2}= zeros(nr2,setup.nvar);
setup.Ss{2}(1,1) =1;    % GDP + 
setup.Ss{2}(2,2)=-1;   % Real Wage -
setup.Ss{2}(3,3)=-1;   % Labor Force Participation -


% shock 3 = domestic labor supply
nr3=4;
setup.Ss{3}= zeros(nr3,setup.nvar);
setup.Ss{3}(1,1)=1;     % GDP + 
setup.Ss{3}(2,2)=-1;   % Real Wage -
setup.Ss{3}(3,3)=1;    % Labor Force Participation + 
setup.Ss{3}(4,4)=-1;   % Immigrants/Participants proxied by Foreign Born over Population


% shock 4 = immigration
nr4=4;
setup.Ss{4}= zeros(nr4,setup.nvar);
setup.Ss{4}(1,1)=1;     % GDP + 
setup.Ss{4}(2,2)=-1;   % Real Wage -
setup.Ss{4}(3,3)= 1;    % Labor Force Participation + 
setup.Ss{4}(4,4)= 1;   % Immigrants/Participants proxied by Foreign Born over Population

