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

import projection.processing.GeometricHash.NodeND;
import projection.util.PointMatch;
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

		try {
			Prefs.savePreferences();
			prepareRegistration(new File(datadir), new File(outputdir));
			register(nLayers);
		} catch(Exception e) {
		IJ.error(e.getMessage());
			e.printStackTrace();
		}
	}

	public static ArrayList<Point3f> getPoints(FindMaxima fm, SphericalMaxProjection smp, short[][] maxima) {
		short[] m = new short[maxima[0].length];
		for(int l = 0; l < maxima.length; l++) {
			for(int i = 0; i < m.length; i++) {
				int mv = m[i] & 0xffff;
				int lv = maxima[l][i] & 0xffff;
				if(lv > mv)
					m[i] = (short)lv;
			}
		}
		short[] ms = m.clone();

		smp.smooth(ms, 3);

		smp.smooth(m, 50);

		for(int i = 0; i < ms.length; i++) {
			int diff = (ms[i] & 0xffff) - (m[i] & 0xffff);
			m[i] = diff > 0 ? (short)diff : 0;
		}
		// System.arraycopy(m, 0, maxima[0], 0, m.length);

		return fm.findMaxima(m, null, 20);
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

	public File getMatrixDirectory() {
		return matrixDirectory;
	}

	private ArrayList<Point3f> tgtPts;

	public static Matrix4f rotateEulerAt(double a1, double a2, double a3, Point3f c) {
		Matrix4f r = new Matrix4f();
		double c1 = Math.cos(a1), s1 = Math.sin(a1);
		double c2 = Math.cos(a2), s2 = Math.sin(a2);
		double c3 = Math.cos(a3), s3 = Math.sin(a3);

		r.m00 = (float)(c3 * c1 - c2 * s1 * s3);
		r.m01 = (float)(-s3 * c1 - c2 * s1 * c3);
		r.m02 = (float)(s2 * s1);
		r.m03 = 0;
		r.m10 = (float)(c3 * s1 + c2 * c1 * s3);
		r.m11 = (float)(-s3 * s1 + c2 * c1 * c3);
		r.m12 = (float)(-s2 * c1);
		r.m13 = 0;
		r.m20 = (float)(s2 * s3);
		r.m21 = (float)(s2 * c3);
		r.m22 = (float)(c2);
		r.m23 = 0;

		Point3f t = new Point3f(c);
		r.transform(c, t);
		r.m03 = c.x - t.x;
		r.m13 = c.y - t.y;
		r.m23 = c.z - t.z;

		r.m30 = r.m31 = r.m32 = 0;
		r.m33 = 1;

		return r;
	}

	double msePrev = -1;

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
			boolean useICP3D = true;
			boolean useICP9D = false;
			if(useICP3D) {
				double mse = ICPRegistration.register(tgtPts, srcPts, mat, smp.center);

				System.out.println("msePrev = " + msePrev + " mse = " + mse);
				double deg = 0;
				double bestMSE = mse;
				Matrix4f bestMat = mat;
//				while(msePrev > 0 && bestMSE > 2 * msePrev && deg < 40) {
//					System.out.println("MSE too bad, try different orientations");
//					// mse got drastically worse, try some other orientations
//
//					deg += 3;
//					double rad = Math.PI * deg / 180;
//					for(int k = -1; k <= 1; k++) {
//						for(int l = -1; l <= 1; l++) {
//							for(int m = -1; m <= 1; m++) {
//								Matrix4f mtmp = rotateEulerAt(k * rad, l * rad, m * rad, smp.center);
//								for(int i = 0; i < srcPts.size(); i++)
//									srcPts.get(i).set(nextTgtPts.get(i));
//								mse = ICPRegistration.register(tgtPts, srcPts, mtmp, smp.center);
//								if(mse < bestMSE) {
//									bestMSE = mse;
//									bestMat = mtmp;
//								}
//							}
//						}
//					}
//				}
				deg = 10;
				double rad = Math.PI * deg / 180;
				for(int k = -1; k <= 1; k++) {
					for(int l = -1; l <= 1; l++) {
						for(int m = -1; m <= 1; m++) {
							Matrix4f mtmp = rotateEulerAt(k * rad, l * rad, m * rad, smp.center);
							for(int i = 0; i < srcPts.size(); i++)
								srcPts.get(i).set(nextTgtPts.get(i));
							mse = ICPRegistration.register(tgtPts, srcPts, mtmp, smp.center);
							if(mse < bestMSE) {
								bestMSE = mse;
								bestMat = mtmp;
							}
						}
					}
				}
				mat = bestMat;
				mse = bestMSE;
				System.out.println("mat = " + mat);
				msePrev = mse;
				System.out.println("msePrev set to " + mse);
			} else if(useICP9D) {
				ArrayList<NodeND> params1 = GeometricHash.hash(tgtPts, false);
				ArrayList<NodeND> params2 = GeometricHash.hash(srcPts, false);
				double mse = ICPRegistrationND.register(params1, params2, mat, smp.center);
				System.out.println("msePrev = " + msePrev + " mse = " + mse);
//				double deg = 0;
//				double bestMSE = mse;
//				Matrix4f bestMat = mat;
//				while(bestMSE > 200 && deg < 40) {
//					System.out.println("MSE too bad, try different orientations");
//					// mse got drastically worse, try some other orientations
//
//					deg += 3;
//					double rad = Math.PI * deg / 180;
//					for(int k = -1; k <= 1; k++) {
//						for(int l = -1; l <= 1; l++) {
//							for(int m = -1; m <= 1; m++) {
//								Matrix4f mtmp = rotateEulerAt(k * rad, l * rad, m * rad, smp.center);
//								for(int i = 0; i < srcPts.size(); i++)
//									srcPts.get(i).set(nextTgtPts.get(i));
//								mse = ICPRegistration.registerND(tgtPts, srcPts, mtmp, smp.center);
//								if(mse < bestMSE) {
//									bestMSE = mse;
//									bestMat = mtmp;
//								}
//							}
//						}
//					}
//				}
//				mat = bestMat;
//				mse = bestMSE;
//				System.out.println("mat = " + mat);
//				msePrev = mse;
//				System.out.println("msePrev set to " + mse);
			} else {
				System.out.println("Found " + tgtPts.size() + " and " + srcPts.size() + " maxima");
				ArrayList<PointMatch> candidates = GeometricHash.establishCorrespondences(tgtPts, srcPts);
				System.out.println(candidates.size() + " candidates left after geometric hashing");
				Matrix4f bestRigid = new Ransac(candidates, smp.center).run();
				if(bestRigid != null)
					mat.set(bestRigid);
				else
					mat.setIdentity();
			}
			overall.mul(mat);
			System.out.println("overall = " + overall);
			TransformIO.saveTransformation(matFile, overall);

		} else {
			System.out.println("loading " + matFile.getAbsolutePath());
			overall = TransformIO.loadTransformation(matFile);
			System.out.println("done");
			msePrev = -1;
		}

		for(int l = 0; l < nLayers; l++) {
			File outputfile = new File(outputDirectory, String.format("%s_%02d.vertices", basename, l));
			if(!outputfile.exists()) {
				System.out.println("applying transformation");
				maxima[l] = smp.applyTransform(overall, maxima[l]);
				System.out.println("saving");
				SphericalMaxProjection.saveShortData(maxima[l], outputfile.getAbsolutePath());
				System.out.println("done");
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
