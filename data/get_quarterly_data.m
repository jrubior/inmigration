%% This file loads the variables that form the VAR model (22 to 26, depending on reduced_form)
function [data,model] = get_quarterly_data(model)

first_obs=model.first_obs;
last_obs = model.last_obs;


vintage_date = model.vintage_date;




%% FRED Data

% Define your FRED API key (replace 'your_api_key' with your actual key)
apiKey = model.apikey;


% Define the ALFRED API URL and series ID
baseUrl   = 'https://api.stlouisfed.org/fred/series/observations';
alfredUrl = 'https://api.alfred.stlouisfed.org/fred/series/observations';


% Define series of interest as either monthly or quarterly, based on the
% level in FRED
% Note: Make sure to use the correct series ID
quarterly_seriesID = {'GDPC1'};



% Initialize a structure to store data
data = table();



% Loop through each monthly series ID
for i = 1:length(quarterly_seriesID)
    % Construct the full API request URL with the frequency parameter
    url = sprintf('%s?series_id=%s&realtime_start=%s&realtime_end=%s&api_key=%s&file_type=json', ...
        baseUrl, quarterly_seriesID{i},vintage_date, vintage_date, apiKey);


    % Fetch the data
    response = webread(url);


    %initialize a table from the downloaded FRED series, remove
    %unneccessary values, convert the series from a string to a number,
    %rename the column according to the series, and clean up the table so
    %you are left with a column of dates and a numeric column labeled
    %according to the series
    df = struct2table(response.observations);

    df.Date = datetime(df.date,'InputFormat','yyyy-MM-dd');
    df = removevars(df, {'date'});

    df = removevars(df, {'realtime_start', 'realtime_end'});
    df.val = str2double(df.value);
    df = renamevars(df,[ 'val'], [sprintf(quarterly_seriesID{i})]);
    df = removevars(df, {'value'});


    %merge the current series into 'data' which is the table that contains
    %all data. Also handles the condition if monthly_seriesID is empty
    if ~isempty(data)
        data = outerjoin(data,df,'MergeKeys',true);
    else
        data = df;

    end

end
%filter dates according to the start and end obs
idx = data.Date<=last_obs & data.Date >= first_obs;
data  = data(idx, :);
model.dates = data.Date;


end