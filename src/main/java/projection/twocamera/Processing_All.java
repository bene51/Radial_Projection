package projection.twocamera;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.Macro;
import ij.plugin.PlugIn;

import java.io.File;
import java.io.IOException;

import javax.vecmath.Matrix4f;

import projection.processing.Map_Projection;
import projection.processing.MultiviewFusion;
import projection.processing.ResampleSMP;
import projection.processing.TimelapseProcessor;

public class Processing_All extends TimelapseProcessor implements PlugIn {

	private File cam1Folder, cam2Folder;
	private File outputFolder;
	private int timepoints, nLayers, nSamples, nChannels;
	private int nProcessed = 0;

	@Override
	public void run(String arg) {
		String options = Macro.getOptions();
		Macro.setOptions(null);
		boolean parallelExecution;

		try {
			cam1Folder   = new File(Macro.getValue(options, "Cam1_folder", ""));
			cam2Folder   = new File(Macro.getValue(options, "Cam2_folder", ""));
			outputFolder = new File(Macro.getValue(options, "Output_folder", ""));
			nSamples     = Integer.parseInt(Macro.getValue(options, "Number_of_samples", "1"));
			nChannels    = Integer.parseInt(Macro.getValue(options, "Number_of_channels", "1"));
			timepoints   = Integer.parseInt(Macro.getValue(options, "Number_of_timepoints", "1"));
			nLayers      = Integer.parseInt(Macro.getValue(options, "Number_of_layers", "1"));
			parallelExecution = Boolean.parseBoolean(Macro.getValue(options, "Execute_timepoints_in_parallel", "false"));
		} catch(Exception e) {
			e.printStackTrace();
		}

		GenericDialogPlus gd = new GenericDialogPlus("Process all");
		gd.addDirectoryField("Cam1_folder", cam1Folder == null ? "" : cam1Folder.getAbsolutePath());
		gd.addDirectoryField("Cam2_folder", cam2Folder == null ? "" : cam2Folder.getAbsolutePath());
		gd.addDirectoryField("Output_folder", outputFolder == null ? "" : outputFolder.getAbsolutePath());

		gd.addNumericField("Number_of_samples", nSamples, 0);
		gd.addNumericField("Number_of_channels", nChannels, 0);
		gd.addNumericField("Number_of_timepoints", timepoints, 0);
		gd.addNumericField("Number_of_layers", nLayers, 0);

		gd.addCheckbox("Execute_timepoints_in_parallel", false);

		gd.showDialog();
		if(gd.wasCanceled())
			return;

		cam1Folder = new File(gd.getNextString());
		cam2Folder = new File(gd.getNextString());
		outputFolder = new File(gd.getNextString());
		nSamples = (int)gd.getNextNumber();
		nChannels = (int)gd.getNextNumber();
		timepoints = (int)gd.getNextNumber();
		nLayers = (int)gd.getNextNumber();
		parallelExecution = gd.getNextBoolean();

		try {
			prepare();
			IJ.log("Start Processing " + timepoints + " timepoints");
			if(parallelExecution)
				processParallel(timepoints);
			else
				processSequential(timepoints);
			IJ.log("done");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private MultiviewFusion[][] fusion;
	private ResampleSMP[][] resample;
	private Map_Projection[][] maps;

	private synchronized void updateProgress() {
		nProcessed++;
		IJ.showStatus(nProcessed + "/" + timepoints);
	}

	public void prepare() throws IOException {
		fusion = new MultiviewFusion[nSamples][nChannels];
		// Register_[] registration = new Register_[nSamples];
		resample = new ResampleSMP[nSamples][nChannels];
		maps = new Map_Projection[nSamples][nChannels];

		for(int s = 0; s < nSamples; s++) {
			for(int c = 0; c < nChannels; c++) {

				File cam1SampleFolder = new File(cam1Folder, "sample" + s);
				File cam2SampleFolder = new File(cam2Folder, "sample" + s);

				IJ.log("Preparing fusion");

				cam1SampleFolder = new File(cam1SampleFolder, "channel" + c);
				cam2SampleFolder = new File(cam2SampleFolder, "channel" + c);

				File fusedFolder = new File(outputFolder, "sample" + s);
				fusedFolder = new File(fusedFolder, "channel" + c);
				if(!fusedFolder.exists())
					fusedFolder.mkdirs();
				fusion[s][c] = new MultiviewFusion();
				fusion[s][c].prepareFusion(new File[] {cam1SampleFolder, cam2SampleFolder}, fusedFolder, nLayers, true);
				fusion[s][c].indicateCameraContributions();
				fusion[s][c].testCameraFusion();

//				IJ.log("Preparing registration");
//				File registeredFolder = new File(fusedFolder, "registered");
//				registration[s] = new Register_();
//				registration[s].prepareRegistration(fusedFolder, registeredFolder);
				File registeredFolder = fusedFolder;

				IJ.log("Preparing resampling");
				File resampledFolder = new File(registeredFolder, "resampled");
				resample[s][c] = new ResampleSMP();
				resample[s][c].prepareResampling(registeredFolder, resampledFolder, nLayers, 2);


				IJ.log("Preparing map projection");
				Matrix4f initial = new Matrix4f();
				initial.setIdentity();
				File mapFolder = new File(registeredFolder, "Mercator");
				maps[s][c] = new Map_Projection();
				maps[s][c].prepareForProjection(nLayers, registeredFolder, initial, Map_Projection.MERCATOR, mapFolder, false, false, false, 1800);
			}
		}
	}

	@Override
	public void processTimepoint(int tp) throws IOException {
		for(int s = 0; s < nSamples; s++) {
			for(int c = 0; c < nChannels; c++) {
				System.out.println("tp " + tp + " sample " + s + " channel " + c);
				long t1 = System.currentTimeMillis();
				fusion[s][c].processTimepoint(tp);
				long t2 = System.currentTimeMillis();
				System.out.println("fused: " + (t2 - t1) + " ms");
//					registration[s].registerTimepoint(tp, nLayers);
				resample[s][c].processTimepoint(tp);
				long t3 = System.currentTimeMillis();
				System.out.println("resampled: " + (t3 - t2) + " ms");
				maps[s][c].processTimepoint(tp);
				long t4 = System.currentTimeMillis();
				System.out.println("projected: " + (t4 - t3) + " ms");
			}
		}
		updateProgress();
	}
}
