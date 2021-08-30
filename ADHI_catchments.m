load_config('config.cfg') 

daily_discharge = readtable([raw_dir filesep 'ADHI_restricted_raw_data' filesep 'daily_discharge.txt']);
coord = readtable([raw_dir filesep 'ADHI_restricted_raw_data' filesep 'coord.txt']);
dates = readtable([raw_dir filesep 'ADHI_restricted_raw_data' filesep 'dates.txt']);
ADHI_code = readtable([raw_dir filesep 'ADHI_restricted_raw_data' filesep 'ADHI_code.txt'],'ReadVariableNames',0,'Delimiter',',');
ADHI_stations = readtable([raw_dir filesep 'ADHI_restricted_raw_data' filesep 'ADHI_stations.csv'],'Delimiter',',');

for ii = 1:size(ADHI_stations,1)
	tic
	
	ADHI_co = ADHI_code.Var1{ii};
	ADHI_co = strrep(ADHI_co,'''','');
	ADHI_co = ADHI_co(6:end);
	
	ID = ADHI_stations.Station_co{ii};
	
	fprintf(['Processing: ' num2str(ii) '    ' num2str(ID) '\n'])
	
	clear DISCHARGE
	DISCHARGE.StationCoords.Lat = coord.Var2(ii);
	DISCHARGE.StationCoords.Lon = coord.Var1(ii);
	DISCHARGE.Station = ADHI_stations.Name{ii};
	refdates = datenum(1900,1,1):datenum(2039,12,31);
	DISCHARGE.Discharge = nan(length(refdates),1,'single');
	adhi_dates = datenum(dates.Var1(:),dates.Var2(:),dates.Var3(:));
	ind = ismember(refdates,adhi_dates);
	DISCHARGE.Discharge(ind) = daily_discharge{:,['Var' num2str(ii)]};
	DISCHARGE.Discharge(DISCHARGE.Discharge<0) = NaN;
	if ~isempty(ADHI_stations.Comment{ii})
		DISCHARGE.Comment = ADHI_stations.Comment{ii};
	end

	clear BOUNDARIES	
	shp_filepath = [raw_dir filesep 'ADHI_restricted_raw_data' filesep 'CatchmentBoundaries' filesep 'ADHIcatch_' ADHI_co '.shp'];
	if ~exist(shp_filepath), continue; end
	shp = shaperead(shp_filepath);
	BOUNDARIES.CatchBounds.Lat = single(shp.Y);
	BOUNDARIES.CatchBounds.Lon = single(shp.X);
	BOUNDARIES.Area = shp.Area;
	
	mkdir([database_dir filesep 'ADHI_restricted_' ID])
	save([database_dir filesep 'ADHI_restricted_' ID filesep 'DISCHARGE.mat'],'DISCHARGE','-v7.3')
	save([database_dir filesep 'ADHI_restricted_' ID filesep 'BOUNDARIES.mat'],'BOUNDARIES','-v7.3')
		
	toc
end