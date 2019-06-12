// DM script to save metadata for the zeta-factor method

void saveMetadata() {
	string data = "{\n";
	data += "\t\"datetime\": \"" + FormatTimeString( GetCurrentTime(),1+16*2 ) + "\",\n";
	data += "\t\"tilt_x\": \"" + EMGetStageAlpha() + " deg\",\n";
	data += "\t\"tilt_y\": \"" + EMGetStageBeta() + " deg\",\n";
	data += "\t\"spot_size\": \"" + (EMGetSpotSize()+1) + "\"\n";
	data += "}";

	string path = "";
	
	Result("Saved metadata: \n" + data + "\n");
	
	SaveAsDialog( "Choose where to save the metadata", "Spectrum x.json", path);
	if (path != "") {
		CreateFile( path )
		number file = OpenFileForWriting( path );
		WriteFile( file, data)
		CloseFile(file);
		result("File saved at: " + path + "\n");
	}
}

saveMetadata();
