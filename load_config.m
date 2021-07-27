function load_config(config_file)

	fid = fopen(config_file,'rt');
	while true
		thisline = fgetl(fid);
		if ~ischar(thisline); break; end  %end of file
		
		thisline = strrep(thisline,' ','');
		thisline = strrep(thisline,'''','');
		
		equal_symbol = strfind(thisline,'=');		
		if length(equal_symbol)==0; continue; end
		
		before = thisline(1:equal_symbol-1);
		after = thisline(equal_symbol+1:end);
		assignin('base',before,after);
	end
end

