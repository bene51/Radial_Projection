package projection.processing;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.PlugIn;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.io.File;
import java.io.FilenameFilter;

public class Create_Overlay extends TimelapseProcessor implements PlugIn {

	public static final int LINECOLOR = new Color(200, 20, 20).getRGB();

	@Override
	public void run(String arg0) {
		GenericDialogPlus gd = new GenericDialogPlus("Create map overlay");
		gd.addDirectoryField("Data directory", "");
		gd.addNumericField("min", 120, 0);
		gd.addNumericField("max", 450, 0);
		gd.addNumericField("nLayers", 1, 0);
		gd.showDialog();
		if(gd.wasCanceled())
			return;

		try {
			File datadir = new File(gd.getNextString());
			if(!datadir.isDirectory()) {
				IJ.error(datadir + " is not a directory");
				return;
			}
			createOverlay(datadir, (int)gd.getNextNumber(), (int)gd.getNextNumber(), (int)gd.getNextNumber());
		} catch(Exception e) {
			e.printStackTrace();
			IJ.error(e.getMessage());
		}


	}

	private File datadir, contribDir, linesDir, outdir;
	private int min, max;
	private int nLayers;

	public void prepareCreateOverlay(File datadir, int min, int max, int nLayers) {
		this.datadir = datadir;
		contribDir = new File(datadir, "contributions");
		linesDir = new File(datadir, "lines");
		outdir = new File(datadir, "overlay");
		outdir.mkdir();
		this.min = min;
		this.max = max;
		this.nLayers = nLayers;
	}


	@Override
	public void processTimepoint(int tp) {
		ImageProcessor contribs = null, li = null;
		String contribFile = new File(contribDir, String.format("tp%04d.tif", tp)).getAbsolutePath();
		try {
			contribs = IJ.openImage(contribFile).getProcessor();
		} catch(Exception e) {
			System.err.println("Could not open contributions file: " + contribFile);
		}

		String linesFile = new File(linesDir, String.format("tp%04d.tif", tp)).getAbsolutePath();
		try {
			li = IJ.openImage(linesFile).getProcessor();
		} catch (Exception e) {
			System.err.println("Could not open lines file");
		}

		for(int l = 0; l < nLayers; l++) {
			File f = new File(datadir, String.format("tp%04d_%02d.tif", tp, l));
			ImageProcessor ip = IJ.openImage(f.getAbsolutePath()).getProcessor();

			ip.setMinAndMax(min, max);
			ip = ip.convertToByte(true);

			ImageProcessor co = contribs != null ? contribs.duplicate() :
				new ColorProcessor(ip.getWidth(), ip.getHeight());

			int len = ip.getWidth() * ip.getHeight();
			for(int i = 0; i < len; i++) {
				int v = ip.get(i);

				// v == 0 : 0.5
				// v == 255: 0
				double c = 1 - (255/2.0 + v/2.0) / 255.0;
				int rgb = co.get(i);
				co.set(i, blend(c, rgb, v));

				if(li != null && li.get(i) == 255)
					co.set(i, LINECOLOR);
			}
			IJ.save(new ImagePlus("", co), new File(outdir, String.format("tp%04d_%02d.tif", tp, l)).getAbsolutePath());
		}
	}


	public static void createOverlay(final File datadir, final int min, final int max, final int nLayers) {
		final Create_Overlay co = new Create_Overlay();
		co.prepareCreateOverlay(datadir, min, max, nLayers);
		final String[] files = datadir.list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.startsWith("tp") && name.endsWith("_00.tif");
			}

		});
		co.processParallel(files.length);
	}

	private static int blend(double f, int rgb, int v) {
		int r = (rgb & 0xff0000) >> 16;
		int g = (rgb & 0xff00) >> 8;
		int b = (rgb & 0xff);
		r = Math.min(255 , (int)Math.round((f * r) + (1-f) * v));
		g = Math.min(255 , (int)Math.round((f * g) + (1-f) * v));
		b = Math.min(255 , (int)Math.round((f * b) + (1-f) * v));
		rgb = (r << 16) + (g << 8) + b;
		return rgb;
	}
}
