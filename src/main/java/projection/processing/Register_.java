package projection.processing;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.Prefs;
import ij.plugin.PlugIn;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import javax.vecmath.Matrix4f;
import javax.vecmath.Point3f;

import projection.util.TransformIO;



public class Register_ implements PlugIn {

	@Override
	public void run(String arg) {
		String datadir = Prefs.get("register_sphere_proj.datadir", "");
		String outputdir = Prefs.get("register_sphere_proj.outputdir", "");
		GenericDialogPlus gd = new GenericDialogPlus("Register sphere projections");
		gd.addDirectoryField("Data directory", datadir);
		gd.addDirectoryField("Output directory", outputdir);
		gd.addNumericField("#layers", 1, 0);
		gd.showDialog();
		if(gd.wasCanceled())
			return;
		datadir = gd.getNextString();
		outputdir = gd.getNextString();
		int nLayers = (int)gd.getNextNumber();

		Prefs.set("register_sphere_proj.datadir", datadir);
		Prefs.set("register_sphere_proj.outputdir", outputdir);
		Prefs.savePreferences();

		try {
			prepareRegistration(new File(datadir), new File(outputdir));
			register(nLayers);
		} catch(Exception e) {
		IJ.error(e.getMessage());
			e.printStackTrace();
		}
	}

	private static ArrayList<Point3f> getPoints(FindMaxima fm, SphericalMaxProjection smp, short[][] maxima) {
		short[] m = new short[maxima[0].length];
		for(int l = 0; l < maxima.length; l++) {
			for(int i = 0; i < m.length; i++) {
				int mv = m[i] & 0xffff;
				int lv = maxima[l][i] & 0xffff;
				if(lv > mv)
					m[i] = (short)lv;
			}
		}
		smp.smooth(m);
		smp.smooth(m);

		return fm.findMaxima(m, 10);
	}

	private SphericalMaxProjection smp;
	private int nVertices;
	private FindMaxima findMaxima;
	private Matrix4f overall;
	private File outputDirectory, dataDirectory, contributionsDirectory, matrixDirectory;
	private int[] contributions;

	public void prepareRegistration(File dDir, File oDir) throws IOException {
		this.dataDirectory = dDir;
		this.outputDirectory = oDir;
		if(!outputDirectory.exists())
			outputDirectory.mkdir();

		File objfile = new File(dataDirectory, "Sphere.obj");
		if(!objfile.exists())
			throw new IllegalArgumentException("Cannot find " + objfile.getAbsolutePath());

		this.contributionsDirectory = new File(outputDirectory, "contributions");
		contributionsDirectory.mkdir();
		this.matrixDirectory = new File(outputDirectory, "transformations");
		matrixDirectory.mkdir();

		this.smp = new SphericalMaxProjection(objfile.getAbsolutePath());
		this.nVertices = smp.getSphere().nVertices;
		smp.saveSphere(new File(outputDirectory, "Sphere.obj").getAbsolutePath());


		contributions = SphericalMaxProjection.loadIntData(new File(dataDirectory, "contributions.vertices").getAbsolutePath(), nVertices);

		this.findMaxima = new FindMaxima(smp);

		this.overall = new Matrix4f();
		overall.setIdentity();
	}

	private ArrayList<Point3f> tgtPts;

	public void registerTimepoint(int tp, int nLayers) throws IOException {
		String basename = String.format("tp%04d", tp);
		File contributionsfile = new File(contributionsDirectory, basename + ".vertices");
		System.out.println("reg: src = " + new File(dataDirectory, basename + ".vertices").getAbsolutePath());
		short[][] maxima = new short[nLayers][];
		for(int l = 0; l < nLayers; l++)
			maxima[l] = SphericalMaxProjection.loadShortData(new File(dataDirectory, String.format("%s_%02d.vertices", basename, l)).getAbsolutePath(), nVertices);

		if(tp == 0) {
			for(int l = 0; l < nLayers; l++) {
				File outputfile = new File(outputDirectory, String.format("%s_%02d.vertices", basename, l));
				SphericalMaxProjection.saveShortData(maxima[l], outputfile.getAbsolutePath());
			}
			SphericalMaxProjection.saveIntData(contributions, contributionsfile.getAbsolutePath());
			tgtPts = getPoints(findMaxima, smp, maxima);
			return;
		}

		ArrayList<Point3f> nextTgtPts = null, srcPts = null;
		String matName = basename + ".matrix";
		File matFile = new File(matrixDirectory, matName);
		if(!matFile.exists()) {
			srcPts = getPoints(findMaxima, smp, maxima);
			if(tgtPts == null) {
				short[][] pmaxima = new short[nLayers][];
				for(int l = 0; l < nLayers; l++)
					pmaxima[l] = SphericalMaxProjection.loadShortData(new File(dataDirectory, String.format("tp%04d_%02d", tp-1, l) + ".vertices").getAbsolutePath(), nVertices);
				tgtPts = getPoints(findMaxima, smp, maxima);
			}

			// make a deep copy of src points, to be used as target points for the next iteration
			nextTgtPts = new ArrayList<Point3f>(srcPts.size());
			for(Point3f p : srcPts)
				nextTgtPts.add(new Point3f(p));

			Matrix4f mat = new Matrix4f();
			mat.setIdentity();
			ICPRegistration.register(tgtPts, srcPts, mat, smp.center);
			overall.mul(mat);
			TransformIO.saveTransformation(matFile, overall);
		} else {
			overall = TransformIO.loadTransformation(matFile);
		}

		for(int l = 0; l < nLayers; l++) {
			File outputfile = new File(outputDirectory, String.format("%s_%02d.vertices", basename, l));
			if(!outputfile.exists()) {
				maxima[l] = smp.applyTransform(overall, maxima[l]);
				SphericalMaxProjection.saveShortData(maxima[l], outputfile.getAbsolutePath());
			}
		}

		if(!contributionsfile.exists()) {
			int[] con = smp.applyTransformNearestNeighbor(overall, contributions);
			SphericalMaxProjection.saveIntData(con, contributionsfile.getAbsolutePath());
		}

		tgtPts = nextTgtPts;
	}

	public void register(int nLayers) throws IOException {
		// obtain list of local maxima files
		HashSet<Integer> tmp = new HashSet<Integer>();
		for(String f : dataDirectory.list())
			if(f.startsWith("tp") && f.endsWith(".vertices"))
				tmp.add(Integer.parseInt(f.substring(2, 6)));
		Integer[] tps = new Integer[tmp.size()];
		tmp.toArray(tps);
		Arrays.sort(tps);

		// register
		for(int i = 0; i < tps.length; i++) {
			registerTimepoint(tps[i], nLayers);
			IJ.showProgress(i, tps.length);
		}
		IJ.showProgress(1);
	}
}
