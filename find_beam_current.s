// ********************************************************************************
// DigitalMicrograph script to find beam current from the CCD camera
// ********************************************************************************
//
// The script is useful when using the zeta-factor method
//
// Author: Daniel Lundeby
// May 2019
//


// Define an wrapper class around Number
Class NumClass : Object
{
	Number var;
	Object init(object self, number v) {
		var = v
		return self;
	}
	Number get(object self) {
		return var;
	}
}


// Calculate the maximum intensity of an image
number calcMax(Image img) {
	number min, max;
	ImageCalculateMinMax( img, 1,  0, min, max );
	return max;
}


// Create a metadata tag as a number
void createTagNumber(TagGroup imgtags, String tagfolder, String tagname, Number value) {
	TagGroup innerTG = TagGroupGetOrCreateTagGroup( imgtags, tagfolder);
	if (! TagGroupDoesTagExist( innerTG, tagname )) {
		TagGroupCreateNewLabeledTag(innerTG, tagname);
	}
	TagGroupSetTagAsNumber(innerTG, tagname, value );
}


// Create a metadata tag as a string
void createTagString(TagGroup imgtags, String tagfolder, String tagname, String value) {
	TagGroup innerTG = TagGroupGetOrCreateTagGroup( imgtags, tagfolder);
	if (! TagGroupDoesTagExist( innerTG, tagname )) {
		TagGroupCreateNewLabeledTag(innerTG, tagname);
	}
	TagGroupSetTagAsString(innerTG, tagname, value );
}


// Create acquisition parameters to be used by the CCD camera
Object createAckParams(Object camera, Number exposure) {
	Object acq_params = CM_GetCameraAcquisitionParameterSet_HighQualityImagingAcquire( camera );
	CM_SetExposure( acq_params, exposure );
	CM_SetProcessing( acq_params, 3 ); // Gain and dark correction
	CM_SetCorrections( acq_params, 1911, 887); // Additional, default corrections
	if ( CM_IsValid_AcquisitionParameters( camera, acq_params ) ) {
		return acq_params;
	}
	else {
		result("Not valid acquisition parameters\n");
	}
}


// Capture a CCD image
Image captureImage(Number exposure) {
	Object camera = CM_GetCurrentCamera();
	Object acq_params = createAckParams(camera, exposure)
	if (! acq_params) {
		return null;
	}
	return CM_AcquireImage( camera, acq_params );
}


// Find the exposure time which causes the max intensity to reach a certain value, by iteration
Number approachMaxIntensity() {
	
	Object aims = alloc(ObjectList);
	AddObjectToList(aims, alloc(NumClass).Init(12000)); // First aim of intensities
	AddObjectToList(aims, alloc(NumClass).Init(12000)); // Second aim of intensities
	AddObjectToList(aims, alloc(NumClass).Init(12000)); // Third aim of intensities
	AddObjectToList(aims, alloc(NumClass).Init(12000)); // Fourth aim of intensities
	
	Number exp = 0.001; //  Initial exposure time, in seconds

	foreach (Object aim; aims) {
		Image img = captureImage(exp);
		result("Exposure time: " + exp + " s, max intensity: " + calcMax(img) + "\n");
		Number factor = aim.get() / calcMax(img);
		exp = exp*factor;
		if (exp > 5) { // Set maximum exposure time in seconds
			exp = 5;
		}
	}
	return exp;

}


// Save metadata to the .dm3 image, and print to console
void saveMetadataAndPrint(Image img, TagGroup imgtags, number index) {
	number counts = sum(img);
	number ma = calcMax(img);
	string datetime = FormatTimeString( GetCurrentTime(), 1+16*2);
	createTagNumber(imgtags, "Acquisition_series:Sum", "Sum_" + index, counts);
	createTagString(imgtags, "Acquisition_series:Time", "Time_" + index, datetime);
	createTagNumber(imgtags, "Acquisition_series:Max", "Max_" + index, ma);
	result("Max intensity: " + ma + ", sum: " + counts + "\n");
}


// Acquire a given number of CCD images
void acquisitionSeries(Number exp, TagGroup imgtags, Number repetitions) {

	result("Acquiring " + repetitions + " images" + "\n");

	Number i = 0;
	while (i < repetitions) {
		Image img = captureImage(exp);
		saveMetadataAndPrint(img, imgtags, i+2)
		i++;
	}

}


// Acquire CCD images as often as possible, for a given number of minutes
void acquisitionSeriesForDuration(Number exp, TagGroup imgtags, Number minutes) {
	
	Number startTime = GetCurrentTime();
	Number endTime = startTime + minutes*60*10000000;
	
	result("Ending at " + FormatTimeString( endTime, 1+16*2) + "\n");
	
	Number i = 0;
	while (GetCurrentTime() < endTime) {
		Image img = captureImage(exp);
		saveMetadataAndPrint(img, imgtags, i+2)
		i++;
	}
}


void run() {
	
	Object camera = CM_GetCurrentCamera();
	if (!CM_GetCameraInserted(camera)) {
		//CM_SetCameraInserted(camera, 1);
		OkDialog("You need to insert the camera before measuring the current\n");
		return;
	}
	
	number cont = ContinueCancelDialog("Confirm that the beam is focused at a hole, and covers most of the CCD.\n\nPress Continue to perform the calculation.");
	if (cont==0) {
		return;
	}
	
	Result("\nStarting exposure series \n");
	Number exp = approachMaxIntensity();
	Result("Ending exposure series. Final exposure: " + exp + "s\n");
	
	Result("Starting acquisition series\n");
	// Capture image for saving, START
	Object acq_params = createAckParams(camera, exp)
	if (! acq_params) {
		return;
	}
	Image img = CM_CreateImageForAcquire( camera, acq_params, "image");
	CM_AcquireImage( camera, acq_params, img );
	TagGroup imgtags = img.imagegettaggroup();
	saveMetadataAndPrint(img, imgtags, 1);
	// Capture image for saving, END
	
	
	// Choose one of the following options:
	acquisitionSeries(exp, imgtags, 9); // Acquire the given number of CCD images
	//acquisitionSeriesForDuration(exp, imgtags, 30); // Acquire CCD images for the given number of minutes
	
	Number saveAfter = 0; // Whether to save a picture after ending acquisition series
	Image imgAfter;
	if (saveAfter) {
		imgAfter = CM_CreateImageForAcquire( camera, acq_params, "image");
		CM_AcquireImage( camera, acq_params, imgAfter );
	}
	Result("Ending acquisition series\n");

	Number current_estimate = Sum(img) / exp / 1.5 * 1.6e-19 * 1e9
	Result("Rough beam current estimate: " + current_estimate + " nA\n");
	
	string path = "";
	SaveAsDialog( "Choose where to save the picture", "current", path);
	if (path != "") {
		SaveAsGatan( img, path );
		result("File saved at: " + path + "\n");
		if (saveAfter) {
			SaveAsGatan( imgAfter, path+"_after" );
			result("File saved at: " + path+"_after" + "\n");
		}
	}
	
}

run();
