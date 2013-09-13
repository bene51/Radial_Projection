package projection;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

import java.io.File;

import projection.processing.Configuration;
import projection.processing.MultiviewRadialMaxProjection;

public class File_MaxProjection implements PlugIn {

	private Configuration conf;
	private int nChannels;
	private File outputdir, inputdir;
	private MultiviewRadialMaxProjection mmsmp;

	public static final String format = "tp%04d_c%03d_a%s" + File.separator + "%04d.tif";

	private File getFile(int tp, int channel, int angleIndex, String[] angleNames, int z) {
		return new File(inputdir, String.format(format, tp, channel, angleNames[angleIndex], z));
	}

	private ImageProcessor getImage(int tp, int channel, int angleIndex, String[] angleNames, int z) {
		String path = getFile(tp, channel, angleIndex, angleNames, z).getAbsolutePath();
		ImagePlus imp = IJ.openImage(path);
		if(imp == null) {
			System.out.println("Cannot open " + path);
			return null;
		}
		return imp.getProcessor();
	}

	@Override
	public void run(String arg) {
		GenericDialogPlus gd = new GenericDialogPlus("File_MaxProjection");
		gd.addFileField("Configuration_file", "");
		gd.addNumericField("nChannels", 1, 0);
		gd.addDirectoryField("Input directory", "");
		gd.addDirectoryField("Output directory", "");
		gd.showDialog();
		if(gd.wasCanceled())
			return;

		try {
			conf = new Configuration(new File(gd.getNextString()));
		} catch(Exception e) {
			IJ.handleException(e);
		}
		nChannels = (int)gd.getNextNumber();
		inputdir = new File(gd.getNextString());
		outputdir = new File(gd.getNextString());

		if(!inputdir.exists() || !inputdir.isDirectory())
			throw new RuntimeException(inputdir + " is not a directory");

		if(!outputdir.exists())
			outputdir.mkdir();

		mmsmp = new MultiviewRadialMaxProjection(outputdir.getAbsolutePath(), conf, nChannels);
		startAcq();
	}

	private void startAcq() {
		for(int c = 0; c < nChannels; c++)
			new File(mmsmp.getOutputDirectory(), "channel" + c).mkdir();

		for(int t = 0; t < conf.nTimepoints; t++) {
			for(int a = 0; a < conf.nAngles; a++) {
				for(int c = 0; c < nChannels; c++) {
					for(int z = 0; z < conf.d; z++) {
						ImageProcessor ip = getImage(t, c, a, conf.angleNames, z);
						mmsmp.process((short[])ip.getPixels(), t, a, c, z);
					}
				}
			}
		}
	}
}
