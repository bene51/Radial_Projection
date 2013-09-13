package projection.processing;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.LutLoader;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import ij.process.ImageProcessor;

import java.awt.image.IndexColorModel;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;

public class Multilayer_Color_Code implements PlugIn {

	@Override
	public void run(String arg) {
		File lutdir = new File(System.getProperty("user.dir"), "luts");
		GenericDialogPlus gd = new GenericDialogPlus("Resample");
		gd.addDirectoryField("Input directory:", "");
		gd.addNumericField("#layers", 1, 0);
		gd.addNumericField("min", 0, 0);
		gd.addNumericField("max", 500, 0);
		String[] luts = lutdir.list();
		gd.addChoice("Lookup table", luts, luts[0]);
		gd.showDialog();
		if(gd.wasCanceled())
			return;

		File datadir = new File(gd.getNextString());
		File outputdir = new File(datadir, "colorcode");
		int nLayers = (int)gd.getNextNumber();
		int min = (int)gd.getNextNumber();
		int max = (int)gd.getNextNumber();
		
		

		try {
			IndexColorModel lut = LutLoader.open(new File(lutdir, gd.getNextChoice()).getAbsolutePath());
			byte[] r = new byte[256];
			byte[] g = new byte[256];
			byte[] b = new byte[256];
			lut.getReds(r);
			lut.getGreens(g);
			lut.getBlues(b);
			colorcode(datadir, outputdir, nLayers, min, max, r, g, b);
		} catch(Exception e) {
			e.printStackTrace();
			IJ.error(e.getMessage());
		}
	}

	private File outputdir;
	private File inputdir;

	public void colorcodeTimepoint(int i, int nLayers, int min, int max) throws IOException {

		File outfile = new File(outputdir, String.format("tp%04d.tif", i));
		ImageStack stack = null;
		for(int l = 0; l < nLayers; l++) {
			File infile = new File(inputdir, String.format("tp%04d_%02d.tif", i, l));
			ImagePlus imp = IJ.openImage(infile.getAbsolutePath());
			ImageProcessor ip = imp.getProcessor();
			ip.setMinAndMax(min, max);
			ip = ip.convertToByte(true);
			ip.setColorModel(colormodel[l]);
			if(l == 0)
				stack = new ImageStack(ip.getWidth(), ip.getHeight());
			stack.addSlice("", ip.convertToRGB());
		}
		
		ZProjector proj = new ZProjector(new ImagePlus("", stack));
		proj.setMethod(ZProjector.MAX_METHOD);
		proj.doRGBProjection();
		ImagePlus projected = proj.getProjection();
		projected.show();
		IJ.save(projected, outfile.getAbsolutePath());
	}
	
	IndexColorModel[] colormodel = null;
	
	private void prepareColorcode(int timepoints, int nlayers, byte[] r, byte[] g, byte[] b) {
		colormodel = new IndexColorModel[nlayers];
		
		byte[] newR = new byte[256];
		byte[] newG = new byte[256];
		byte[] newB = new byte[256];

		for (int i = 0; i < nlayers; i++) {
			int colorscale = (int)Math.floor((256.0 / nlayers) * i);
			for (int j = 0; j < 256; j++) {
				double intensityfactor = j / 255.0;
				newR[j] = (byte)Math.round((r[colorscale] & 0xff) * intensityfactor);
				newG[j] = (byte)Math.round((g[colorscale] & 0xff) * intensityfactor);
				newB[j] = (byte)Math.round((b[colorscale] & 0xff) * intensityfactor);
			}
			colormodel[i] = new IndexColorModel(8, 256, newR, newG, newB);
		}
	}

	public static void colorcode(File inputdir, File outputdir, int nLayers, int min, int max, byte[] r, byte[] g, byte[] b) throws IOException {
		Multilayer_Color_Code res = new Multilayer_Color_Code();
		
		res.inputdir = inputdir;
		res.outputdir = outputdir;
		
		if(!res.outputdir.exists())
			res.outputdir.mkdir();

		
		// obtain list of local maxima files
		HashSet<Integer> tmp = new HashSet<Integer>();
		for(String file : inputdir.list())
			if(file.startsWith("tp") && file.endsWith(".tif"))
				tmp.add(Integer.parseInt(file.substring(2, 6)));
		Integer[] tps = new Integer[tmp.size()];
		tmp.toArray(tps);
		Arrays.sort(tps);
		
		res.prepareColorcode(tps.length, nLayers, r, g, b);

		// register
		for(int i = 0; i < tps.length; i++) {
			res.colorcodeTimepoint(tps[i], nLayers, min, max);
			IJ.showProgress(i, tps.length);
		}
		IJ.showProgress(1);
	}
}
