package projection.processing;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.plugin.PlugIn;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;

import org.scijava.vecmath.Point3f;

public class TV_Filter implements PlugIn {

	@Override
	public void run(String arg0) {
		GenericDialogPlus gd = new GenericDialogPlus("TV Filter");
		gd.addDirectoryField("Data directory", "");

		gd.addNumericField("lambda", 0.1, 5);
		gd.addNumericField("tolerance", 1, 5);
		gd.addNumericField("a", 0.0001, 5);
		gd.showDialog();
		if(gd.wasCanceled())
			return;

		String datadir = gd.getNextString();
		double lambda = gd.getNextNumber();
		double tolerance = gd.getNextNumber();
		double a = gd.getNextNumber();

		try {
			process(datadir, lambda, tolerance, a);
		} catch (IOException e) {
			e.printStackTrace();
			IJ.error(e.getMessage());
		}
	}

	public void process(String dir, double lambda, double tolerance, double a) throws IOException {
		File dataf = new File(dir);
		File outdir = new File(dir, "filtered");
		outdir.mkdir();
		String objfile = new File(dataf, "Sphere.obj").getAbsolutePath();
		SphericalMaxProjection smp = new SphericalMaxProjection(objfile);
		Point3f[] vertices = smp.getSphere().getVertices();
		int[] faces = smp.getSphere().getFaces();
		int[][] neighbors = calculateNeighbors(faces, vertices.length);

		File[] files = dataf.listFiles();
		for(int i = 0; i < files.length; i++) {
			if(!files[i].getName().endsWith(".vertices"))
				continue;

			short[] original = smp.loadMaxima(files[i].getAbsolutePath());
			double diff = 0;
			float[] current = new float[original.length];
			for(int c = 0; c < current.length; c++)
				current[c] = original[c];
			int iter = 0;
			do {
				float[] next = new float[vertices.length];
				diff = step(neighbors, current, next, original, lambda, a);
				current = next;
				iter++;
			} while(diff > tolerance);

			for(int c = 0; c < current.length; c++)
				original[c] = (short)Math.round(current[c]);
			System.out.println(iter + " iterations");

			smp.saveMaxima(new File(outdir, files[i].getName()).getAbsolutePath(), original);
		}
	}

	@SuppressWarnings("unchecked")
	private int[][] calculateNeighbors(int[] faces, int nVertices) {
		HashSet<Integer>[] set = new HashSet[nVertices];
		for(int i = 0; i < nVertices; i++)
			set[i] = new HashSet<Integer>();

		for(int i = 0; i < faces.length; i += 3) {
			int f1 = faces[i];
			int f2 = faces[i + 1];
			int f3 = faces[i + 2];
			set[f1].add(f2);
			set[f1].add(f3);

			set[f2].add(f1);
			set[f2].add(f3);

			set[f3].add(f1);
			set[f3].add(f2);
		}

		int[][] neigh = new int[nVertices][];
		for(int i = 0; i < nVertices; i++) {
			neigh[i] = new int[set[i].size()];
			int n = 0;
			for(int neighbor : set[i])
				neigh[i][n++] = neighbor;
		}
		return neigh;
	}

	private double step(int[][] neighbors, float[] u, float[] unext, short[] orig, double lambda, double a) {
		float[] lv = calculateLocalVariation(neighbors, u, a);

		double diff = 0;
		for(int v = 0; v < neighbors.length; v++) {
			double w_ag = 0;
			for(int n : neighbors[v])
				w_ag += 1 / lv[n] + 1 / lv[v];

			double h_aa = lambda / (lambda + w_ag);
			double F = h_aa * (orig[v] & 0xffff);

			for(int n : neighbors[v]) {
				double w_ab = 1.0 / lv[v] + 1.0 / lv[n];
				double h_ab = w_ab / (lambda + w_ag);
				F += h_ab * u[n];
			}
			unext[v] = (float)F;
			diff = Math.max(diff, Math.abs(u[v] - F));
		}
		return diff;
	}

	private float[] calculateLocalVariation(int[][] neighbors, float[] data, double a) {
		float[] newdata = new float[data.length];
		for(int v = 0; v < neighbors.length; v++) {
			double sum = a * a;
			for(int n : neighbors[v])
				sum += Math.pow(data[n] - data[v], 2);
			newdata[v] = (float)Math.sqrt(sum);
		}
		return newdata;
	}
}
