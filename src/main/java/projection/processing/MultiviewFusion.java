package projection.processing;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.plugin.PlugIn;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Set;
import java.util.TreeSet;

import org.scijava.vecmath.Matrix4f;
import org.scijava.vecmath.Point3f;

import projection.util.TransformIO;


public class MultiviewFusion extends TimelapseProcessor implements PlugIn {

	static final int[] colors = new int[] {
		new Color(206, 47, 42).getRGB(),
		new Color(241, 183, 51).getRGB(),
		new Color(162, 198, 231).getRGB(),
		new Color(42, 76, 149).getRGB()
	};


	private InputFolder inputdir;

	private File outputdir;
	private SphericalMaxProjection smp;
	private float[][] weights;

	// the forward transforms, i.e. the transformation that transforms
	// the reference sphere to the individual views.
	private Matrix4f[] transforms;
	private boolean saveOutput;
	private double[] angles;
	private int[][] usedVertexIndices;
	private String[] angleNames;
	private int nLayers;

	@Override
	public void run(String args) {
		GenericDialogPlus gd = new GenericDialogPlus("Multiview fusion");
		gd.addDirectoryField("Input folder", "");
		gd.addNumericField("#layers", 1, 0);
		gd.showDialog();
		if(gd.wasCanceled())
			return;

		String folder = gd.getNextString();
		int nLayers = (int)gd.getNextNumber();

		try {
			fuse(new File(folder), nLayers, true);
		} catch (IOException e) {
			IJ.handleException(e);
		}
	}

	public static void fuse(File folder, int nLayers, boolean saveOutput) throws IOException {

		final MultiviewFusion tcf = new MultiviewFusion();
		tcf.prepareFusion(folder, nLayers, saveOutput);
		tcf.indicateCameraContributions();
		tcf.testCameraFusion();

		final Set<Integer> tps = new TreeSet<Integer>();
		for(File f : folder.listFiles()) {
			String name = f.getName();
			if(name.startsWith("tp") && name.endsWith(".vertices")) {
				int n = Integer.parseInt(name.substring(2, 6));
				tps.add(n);
			}
		}


		int nTimepoints = tps.size();
		tcf.processParallel(nTimepoints);
	}


	/**
	 * Prepare fusion, using a single input folder.
	 * @param saveOutput
	 * @throws IOException
	 */
	public void prepareFusion(File folder, int nLayers, boolean saveOutput) throws IOException {
		this.inputdir = new SingleInputFolder(folder);
		this.saveOutput = saveOutput;
		File objfile = new File(folder, "Sphere.obj");
		// TODO this is probably always the case
		if(!objfile.exists())
			objfile = new File(folder.getParentFile(), "Sphere.obj");
		this.smp = new SphericalMaxProjection(objfile.getAbsolutePath());
		this.outputdir = new File(folder, "fused");
		this.nLayers = nLayers;

		String[] angleFiles = folder.list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(".mat");
			}

		});
		Arrays.sort(angleFiles);
		int nAngles = angleFiles.length;

		this.angles = new double[nAngles];
		this.angleNames = new String[nAngles];
		double[] apertures = new double[nAngles];
		this.transforms = new Matrix4f[nAngles];
		this.usedVertexIndices = new int[nAngles][];


		for(int aIndex = 0; aIndex < nAngles; aIndex++) {
			File f = new File(folder, angleFiles[aIndex]);
			BufferedReader in = new BufferedReader(new FileReader(f));
			angleNames[aIndex] = in.readLine();
			angles[aIndex] = Double.parseDouble(in.readLine());
			// colors[aIndex] = cm.getRGB(aIndex);
			apertures[aIndex] = Double.parseDouble(in.readLine());
			Matrix4f m = TransformIO.fromString(in.readLine());
			transforms[aIndex] = m;
			in.close();

			f = new File(folder, angleFiles[aIndex].replace(".mat", ".indices"));
			usedVertexIndices[aIndex] = MaximumProjector.loadVertexIndices(f);
		}

		if(saveOutput) {
			if(!outputdir.exists())
				outputdir.mkdirs();
			smp.saveSphere(new File(outputdir, "Sphere.obj").getAbsolutePath());
		}

		Point3f center = smp.getCenter();
		weights = new float[angles.length][smp.getSphere().nVertices];
		for(int a = 0; a < angles.length; a++) {
			Point3f cen = new Point3f(center);
			if(transforms[a] != null)
				transforms[a].transform(cen);
			AngleWeighter aw = new AngleWeighter(AngleWeighter.X_AXIS,  angles[a], apertures[a], new Point3f(cen));
			for(int v = 0; v < smp.getSphere().nVertices; v++) {
				Point3f vtx = smp.getSphere().vertices[v];
				if(transforms[a] != null) {
					vtx = new Point3f(vtx);
					transforms[a].transform(vtx);
				}
				weights[a][v] = aw.getWeight(vtx.x, vtx.y, vtx.z);
			}
		}
	}

	// one file for each angle
	public void prepareFusion(File[] folders, File outputdir, int nLayers, boolean saveOutput) throws IOException {
		this.saveOutput = saveOutput;
		File objfile = new File(folders[0], "Sphere.obj");
		// TODO this is probably always the case
		if(!objfile.exists())
			objfile = new File(folders[0].getParentFile(), "Sphere.obj");
		this.smp = new SphericalMaxProjection(objfile.getAbsolutePath());
		this.outputdir = outputdir;
		this.nLayers = nLayers;

		ArrayList<File> angleFiles = new ArrayList<File>();
		for(int i = 0; i < folders.length; i++) {
			angleFiles.addAll(Arrays.asList(folders[i].listFiles(new FilenameFilter() {
				@Override
				public boolean accept(File dir, String name) {
					return name.endsWith(".mat");
				}

			})));
			Collections.sort(angleFiles);
		}
		int nAngles = angleFiles.size();

		this.angles = new double[nAngles];
		this.angleNames = new String[nAngles];
		double[] apertures = new double[nAngles];
		this.transforms = new Matrix4f[nAngles];
		this.usedVertexIndices = new int[nAngles][];


		for(int aIndex = 0; aIndex < nAngles; aIndex++) {
			File f = angleFiles.get(aIndex);
			BufferedReader in = new BufferedReader(new FileReader(f));
			angleNames[aIndex] = in.readLine();
			angles[aIndex] = Double.parseDouble(in.readLine());
			apertures[aIndex] = Double.parseDouble(in.readLine());
			transforms[aIndex] = TransformIO.fromString(in.readLine());
			in.close();
		}

		folders = assignFoldersToAngles(folders);
		this.inputdir = new MultipleInputFolder(folders);

		for(int aIndex = 0; aIndex < nAngles; aIndex++) {
			File f = new File(folders[aIndex], angleFiles.get(aIndex).getName().replace(".mat", ".indices"));
			usedVertexIndices[aIndex] = MaximumProjector.loadVertexIndices(f);
		}

		if(saveOutput) {
			if(!outputdir.exists())
				outputdir.mkdirs();
			smp.saveSphere(new File(outputdir, "Sphere.obj").getAbsolutePath());
		}

		Point3f center = smp.getCenter();
		weights = new float[angles.length][smp.getSphere().nVertices];
		for(int a = 0; a < angles.length; a++) {
			Point3f cen = new Point3f(center);
			if(transforms[a] != null)
				transforms[a].transform(cen);
			AngleWeighter aw = new AngleWeighter(AngleWeighter.X_AXIS,  angles[a], apertures[a], new Point3f(cen));
			for(int v = 0; v < smp.getSphere().nVertices; v++) {
				Point3f vtx = smp.getSphere().vertices[v];
				if(transforms[a] != null) {
					vtx = new Point3f(vtx);
					transforms[a].transform(vtx);
				}
				weights[a][v] = aw.getWeight(vtx.x, vtx.y, vtx.z);
			}
		}
	}

	private File[] assignFoldersToAngles(File[] folders) {
		// tp0000_0045_l00.vertices
		File[] ret = new File[angleNames.length];
		for(int a = 0; a < angleNames.length; a++) {
			String aName = angleNames[a];
			for(File dir : folders) {
				File f = new File(dir, "tp0000_" + aName + "_l00.vertices");
				if(f.exists()) {
					ret[a] = dir;
					break;
				}
			}
		}
		return ret;
	}

	public void testCameraFusion() throws IOException {
		File dir = new File(outputdir, "test");
		dir.mkdir();

		Point3f[] vertices = smp.getSphere().getVertices();
		for(int a = 0; a < angles.length; a++) {
			File outf = new File(dir, getFile(0, a, 0));
			if(outf.exists())
				continue;
			short[] res = new short[vertices.length];
			for(int v = 0; v < vertices.length; v++)
				res[v] = (short)(100 * weights[a][v]);
			SphericalMaxProjection.saveShortData(res, outf.getAbsolutePath());
		}
	}

	public void indicateCameraContributions() throws IOException {
		File out = new File(outputdir, "contributions.vertices");
		if(out.exists())
			return;

		Point3f[] vertices = smp.getSphere().getVertices();
		int[] res = new int[vertices.length];
		for(int v = 0; v < vertices.length; v++) {
			double r = 0, g = 0, b = 0;
			double sum = 0;
			for(int a = 0; a < angles.length; a++) {
				float w = weights[a][v];
				int c = colors[a % colors.length];

				int red   = (c & 0xff0000) >> 16;
				int green = (c & 0xff00) >> 8;
				int blue  = (c & 0xff);

				r += w * red;
				g += w * green;
				b += w * blue;

				sum += w;
			}
			r /= sum;
			g /= sum;
			b /= sum;

			int ir = r > 255 ? 255 : (int)r;
			int ig = g > 255 ? 255 : (int)g;
			int ib = b > 255 ? 255 : (int)b;
			res[v] = (ir << 16) + (ig << 8) + ib;
		}
		SphericalMaxProjection.saveIntData(res, out.getAbsolutePath());
	}

	private String getFile(int tp, int aIndex, int layer) {
		return MultiviewRadialMaxProjection.getFileName(tp, angleNames, aIndex, layer);
	}

	@Override
	public void processTimepoint(int tp) throws IOException {
		boolean done = true;
		for(int l = 0; l < nLayers; l++) {
			if(!new File(outputdir, String.format("tp%04d_%02d.vertices", tp, l)).exists()) {
				done = false;
				break;
			}
		}
		if(done)
			return;

		int nVertices = smp.getSphere().nVertices;

		short[][][] m = new short[angles.length][nLayers][nVertices];
		for(int a = 0; a < angles.length; a++) {
			for(int l = 0; l < nLayers; l++) {
				// m[a][l] = SphericalMaxProjection.loadShortData(inputdir.getFile(a, getFile(tp, a, l)).getAbsolutePath(),  nVertices);
				File f = inputdir.getFile(a, getFile(tp, a, l));
				MaximumProjector.loadVertices(f, usedVertexIndices[a], m[a][l]);
			}
		}

		short[] sData = new short[nVertices];
		Point3f[] vertices = smp.getSphere().getVertices();

		for(int l = 0; l < nLayers; l++) {
			File out = new File(outputdir, String.format("tp%04d_%02d.vertices", tp, l));
			float[] res = new float[vertices.length];
//			for(int v = 0; v < vertices.length; v++) {
//				float sum = 0;
//				for(int a = 0; a < angles.length; a++) {
//					float w = weights[a][v];
//					float ma = (m[a][l][v] & 0xffff);
//					sum += w;
//					res[v] += (w * ma);
//				}
//				res[v] = sum == 0 ? 0 : res[v] / sum;
//			}
			for(int v = 0; v < vertices.length; v++) {
				float max = 0;
				for(int a = 0; a < angles.length; a++) {
					float ma = (m[a][l][v] & 0xffff);
					if(ma > max)
						max = ma;
				}
				res[v] = max;
			}
			for(int v = 0; v < res.length; v++)
				sData[v] = (short)res[v];

			if(saveOutput)
				SphericalMaxProjection.saveShortData(sData, out.getAbsolutePath());
		}
	}

	static interface InputFolder {
		File getFile(int indexInFolderArray, String filename);
	}

	static class SingleInputFolder implements InputFolder {

		private final File folder;

		SingleInputFolder(File folder) {
			this.folder = folder;
		}

		@Override
		public File getFile(int indexInFolderArray, String file) {
			return new File(folder, file);
		}
	}

	static class MultipleInputFolder implements InputFolder {

		private final File[] folders;

		MultipleInputFolder(File[] folders) {
			this.folders = folders;
		}

		@Override
		public File getFile(int indexInFolderArray, String file) {
			File f = new File(folders[indexInFolderArray], file);
			while(!f.exists()) {
				try {
					Thread.sleep(5000);
				} catch(InterruptedException e) {
					e.printStackTrace();
				}
			}
			return f;
		}
	}
}
