%% Housekeeping
clear variables;clc;
addpath('td'); % useful for interpolation


% user specifications
model_m.apikey ='';% your FRED API goes here
model_q.apikey = model_m.apikey;


model_m.vintage_date = '2026-04-01';
model_q.vintage_date = model_m.vintage_date;

model_m.first_obs = datetime(1959, 01, 01);
model_m.last_obs  = datetime(2025, 12, 01);

model_q.first_obs = datetime(1959, 01, 01);
model_q.last_obs  = datetime(2025, 12, 01);

model.first_obs = datetime(1994, 01, 01); % DO NOT CHANGE
model.last_obs  = datetime(2019, 12, 01);


[data_m,model_m] = get_monthly_data(model_m);
[data_q,model_q] = get_quarterly_data(model_q);




%% Interpolate quarterly GDP with industrial production and real personal consumption expenditures
% we use the entire sample since 1959
GDPC1 = data_q.GDPC1;
IP    = data_m.INDPRO;
RPCE  = data_m.DPCERA3M086SBEA;

Y  = GDPC1;     % Y: Nx1 ---> vector of low frequency data
x  = [IP,RPCE]; % x: nxp ---> matrix of high frequency indicators (without intercept)
ta = 2;     % type of disaggregation ---> average (index)
sc = 3;     % quarterly to monthly
type = 0;   % estimation method:  (0) weighted least squares  (1) maximum likelihood
opC  = 1;   % no intercept in hf model

res = chowlin(Y,x,ta,sc,type,opC,[]);

MGDP = res.y;


%% Real Wage
REALWAGE = data_m.AHETPI./data_m.PCEPI;


%% Participation rate
idx_shutdown = find(data_m.Date == datetime(2025,10,1));
data_m.CIVPART(idx_shutdown,1) = (data_m.CIVPART(idx_shutdown-1,1)+data_m.CIVPART(idx_shutdown+1,1))/2;
LFPART = data_m.CIVPART;


%% Foreign born to population ratio
%data_m.LNU00073395(idx_shutdown,1) = (data_m.LNU00073395(idx_shutdown-1,1)+data_m.LNU00073395(idx_shutdown+1,1))/2;
%FBORNPOP = 100*data_m.LNU00073395./data_m.POPTHM;
mandelman_data = readtable('monthly_foreign_born_share.csv');
idx            = find(mandelman_data.YEAR == year(model.last_obs) & mandelman_data.MONTH == month(model.last_obs));
FBORNPOP       = 100*mandelman_data.foreign_born_share(1:idx);



% Unemployment rate
UNRATE = data_m.UNRATE;


%% Dataset
idx_first = find(data_m.Date == model.first_obs);
idx_last = find(data_m.Date == model.last_obs);

dataset = [100*log(MGDP(idx_first:idx_last,1)),100*log(REALWAGE(idx_first:idx_last,1)),LFPART(idx_first:idx_last,1),FBORNPOP,UNRATE(idx_first:idx_last,1)];


writematrix(dataset, 'data_FURLANETTO_ROBSTAD_USA.csv') 