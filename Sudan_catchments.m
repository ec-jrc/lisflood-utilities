load_config('config.cfg') 

Streamflow_Guages = shaperead([raw_dir filesep 'Sudan_restricted_raw_data' filesep 'GIS Layers' filesep 'Streamflow Guages.shp']);
Metadata = readtable([raw_dir filesep 'Sudan_restricted_raw_data' filesep 'Sudan Row Streamflow Data' filesep 'Sudan_Streamflow_Data_(1965-2020).xlsx'],...
	'Sheet','Metadata');
Observed_Daily_Streamflow_MCM = readtable([raw_dir filesep 'Sudan_restricted_raw_data' filesep 'Sudan Row Streamflow Data' filesep 'Sudan_Streamflow_Data_(1965-2020).xlsx'],...
	'Sheet','Observed_Daily_Streamflow_MCM');
	
for ii = 1:size(Metadata,1)
	tic
	
	id = Metadata.ID(ii);
	if isnan(id)
		ID = num2str(ii);
	else	
		ID = num2str(Metadata.ID(ii));
	end
	
	fprintf(['Processing: ' num2str(ii) '    ' num2str(ID) '\n'])
	
	clear DISCHARGE
	DISCHARGE.StationCoords.Lat = Metadata.N(ii);
	DISCHARGE.StationCoords.Lon = Metadata.E(ii);
	DISCHARGE.Station = [Metadata.River{ii} ' at ' Metadata.Station{ii}];
	refdates = datenum(1900,1,1):datenum(2039,12,31);
	DISCHARGE.Discharge = nan(length(refdates),1,'single');
	Sudan_dates = datenum(Observed_Daily_Streamflow_MCM.Date(:));
	ind = ismember(refdates,Sudan_dates);
	station_name = strrep(strrep(strrep(Metadata.Station{ii},' ',''),')','_'),'(','_');	
	qobs = Observed_Daily_Streamflow_MCM{:,[station_name]};
	if iscell(qobs(1))
		DISCHARGE.Discharge(ind) = str2double(qobs);
	else
		DISCHARGE.Discharge(ind) = qobs;
	end
	DISCHARGE.Discharge(DISCHARGE.Discharge<0) = NaN;

	% Remove suspect discharge values
	if strcmp(Metadata.Station{ii},'Menagil (Main Canal)')
		DISCHARGE.Discharge(DISCHARGE.Discharge>50) = NaN;
	elseif strcmp(Metadata.Station{ii},'Girba (Main Canal)')
		DISCHARGE.Discharge(DISCHARGE.Discharge>12) = NaN;
	elseif strcmp(Metadata.Station{ii},'Malakal')
		DISCHARGE.Discharge(DISCHARGE.Discharge>150) = NaN;
	end
	
	%clear BOUNDARIES	
	%shp_filepath = [raw_dir filesep 'Sudan_restricted_raw_data' filesep 'CatchmentBoundaries' filesep 'Sudancatch_' Sudan_co '.shp'];
	%if ~exist(shp_filepath), continue; end
	%shp = shaperead(shp_filepath);
	%BOUNDARIES.CatchBounds.Lat = single(shp.Y);
	%BOUNDARIES.CatchBounds.Lon = single(shp.X);
	%BOUNDARIES.Area = shp.Area;
	
	
	
	mkdir([database_dir filesep 'Sudan_restricted_' ID])
	save([database_dir filesep 'Sudan_restricted_' ID filesep 'DISCHARGE.mat'],'DISCHARGE','-v7.3')
	%save([database_dir filesep 'Sudan_restricted_' ID filesep 'BOUNDARIES.mat'],'BOUNDARIES','-v7.3')
		
	toc
end