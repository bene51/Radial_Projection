package projection.processing;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.plugin.PlugIn;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;

import javax.vecmath.Point3f;

public class ResampleSMP extends TimelapseProcessor implements PlugIn {

	@Override
	public void run(String arg) {
		GenericDialogPlus gd = new GenericDialogPlus("Resample");
		gd.addDirectoryField("Input directory:", "");
		gd.addNumericField("Levels: ", 2, 3);
		gd.addNumericField("#layers", 1, 0);
		gd.showDialog();
		if(gd.wasCanceled())
			return;

		File datadir = new File(gd.getNextString());
		File outputdir = new File(datadir, "resampled");
		int f = (int)gd.getNextNumber();
		int nLayers = (int)gd.getNextNumber();

		try {
			resample(datadir, outputdir, f, nLayers);
		} catch(Exception e) {
			e.printStackTrace();
			IJ.error(e.getMessage());
		}
	}

	private File outputdir;
	private File inputdir;
	private File contribInDir, contribOutDir;
	private SphericalMaxProjection smp;
	private int levels;
	private int nLayers;
	private IndexedTriangleMesh sphere;
	private int[] indices;

	// used if resampling unfused data:
	private int[][] usedVertexIndices = null;
	private String[] angleNames = null;

	public void prepareResampling(File inputdir, File outputdir, int nLayers, int levels) throws IOException {
		this.inputdir = inputdir;
		this.outputdir = outputdir;
		this.levels = levels;
		this.nLayers = nLayers;
		this.smp = new SphericalMaxProjection(new File(inputdir, "Sphere.obj").getAbsolutePath());

		final int n = smp.getSphere().nVertices;
		final float radius = smp.getRadius();
		Point3f center = smp.getCenter();

		final int subd = (int)((-20 + Math.sqrt(400 - 40 * (12 - n))) / 20 + 1);
		final int oldlevels = (int)Math.floor(Math.log(subd) / Math.log(2));
		final int newlevels = oldlevels - levels;

		IJ.showStatus("Creating icosahedron");
		Icosahedron icosa = new Icosahedron(radius);

		IJ.showStatus("Creating buckyball");
		sphere = icosa.createBuckyballFast(radius, newlevels);
		for(Point3f p : sphere.getVertices())
			p.add(center);

		IJ.showStatus("Creating indices");
		indices = new int[sphere.nVertices];
		for(int i = 0; i < indices.length; i++)
			indices[i] = smp.getNearestNeighbor(sphere.getVertices()[i]);

		if(!outputdir.exists())
			outputdir.mkdir();

		contribInDir = new File(inputdir, "contributions");
		contribOutDir = new File(outputdir, "contributions");
		contribOutDir.mkdir();
		SphericalMaxProjection.saveSphere(sphere, new File(outputdir, "Sphere.obj").getAbsolutePath());

		// resample contributions if present
		File contfile = new File(inputdir, "contributions.vertices");
		if(contfile.exists()) {
			int[] overlayold = SphericalMaxProjection.loadIntData(contfile.getAbsolutePath(), smp.getSphere().nVertices);
			int[] overlaynew = new int[sphere.nVertices];
			for(int j = 0; j < indices.length; j++)
				overlaynew[j] = overlayold[indices[j]];
			SphericalMaxProjection.saveIntData(overlaynew, new File(outputdir, "contributions.vertices").getAbsolutePath());
		}

		String[] angleFiles = inputdir.list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(".indices");
			}

		});
		if(angleFiles != null && angleFiles.length > 0) {
			Arrays.sort(angleFiles);
			int nAngles = angleFiles.length;
			this.usedVertexIndices = new int[nAngles][];
			this.angleNames = new String[nAngles];

			for(int aIndex = 0; aIndex < nAngles; aIndex++) {
				File f = new File(inputdir, angleFiles[aIndex]);
				usedVertexIndices[aIndex] = MaximumProjector.loadVertexIndices(f);
				angleNames[aIndex] = angleFiles[aIndex].substring(0, 4);
			}
		}
	}

	@Override
	public void processTimepoint(final int tp) throws IOException {
		int nAngles = usedVertexIndices != null ? usedVertexIndices.length : 1;

		for(int a = 0; a < nAngles; a++) {
			for(int l = 0; l < nLayers; l++) {
				String basename = usedVertexIndices == null ?
						String.format("tp%04d_%02d", tp, l) :
						String.format("tp%04d_%s_l%02d", tp, angleNames[a], l);
				String name = basename + ".vertices";

				File infile = new File(inputdir, name);
				File outfile = new File(outputdir, name);
				if(outfile.exists())
					continue;

				short[] maxima = null;
				if(usedVertexIndices == null) {
					maxima = smp.loadMaxima(infile.getAbsolutePath());
				} else {
					maxima = new short[smp.getSphere().nVertices];
					MaximumProjector.loadVertices(infile, usedVertexIndices[a], maxima);
				}

				for(int j = 0; j < levels; j++)
					smp.smooth(maxima);

				short[] resampled = new short[sphere.nVertices];
				for(int j = 0; j < indices.length; j++)
					resampled[j] = maxima[indices[j]];

				SphericalMaxProjection.saveShortData(resampled, outfile.getAbsolutePath());
			}
		}

		String name = String.format("tp%04d_%02d", tp, 0);
		File contfile = new File(contribInDir, name);
		if(!contfile.exists())
			return;
		int[] overlayold = SphericalMaxProjection.loadIntData(contfile.getAbsolutePath(), smp.getSphere().nVertices);
		int[] overlaynew = new int[sphere.nVertices];
		for(int j = 0; j < indices.length; j++)
			overlaynew[j] = overlayold[indices[j]];
		SphericalMaxProjection.saveIntData(overlaynew, new File(contribOutDir, name).getAbsolutePath());
	}

	/*
	 * Overall no of vertices = 12 * (Sum_{i=1}^{subd + 2}(i) - 3) + 12 - 30 * subd;
	 *                        =
	 */
	public static void resample(File inputdir, File outputdir, int f, int nLayers) throws IOException {
		ResampleSMP res = new ResampleSMP();
		res.prepareResampling(inputdir, outputdir, nLayers, f);

		// obtain list of local maxima files
		HashSet<Integer> tmp = new HashSet<Integer>();
		for(String file : inputdir.list())
			if(file.startsWith("tp") && file.endsWith(".vertices"))
				tmp.add(Integer.parseInt(file.substring(2, 6)));
		int nTimepoints = tmp.size();
		res.processParallel(nTimepoints);
	}
}
