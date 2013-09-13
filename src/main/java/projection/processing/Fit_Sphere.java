package projection.processing;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.StackConverter;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;

import javax.vecmath.Point3f;
import javax.vecmath.Point4f;

import Jama.Matrix;
import customnode.MeshMaker;

public class Fit_Sphere implements PlugInFilter {
	private ImagePlus image;
	private double x0, y0, z0, r;

	public Fit_Sphere() {}

	public Fit_Sphere(ImagePlus image) {
		this.image = image;
	}

	@Override
	public int setup(String arg, ImagePlus image) {
		this.image = image;
		return DOES_8G | DOES_16 | DOES_32;
	}

	public void getCenter(Point3f center) {
		center.x = (float)x0;
		center.y = (float)y0;
		center.z = (float)z0;
	}

	public double getRadius() {
		return r;
	}

	@Override
	public void run(ImageProcessor ip) {
		double threshold = IJ.getNumber(
				"Threshold", 127);
		if(threshold == IJ.CANCELED)
			return;
		fit(threshold);
		IJ.log("x0 = " + x0 + "\ny0 = " + y0 + "\nz0 = " + z0 + "\nr = " + r);
		getControlImage().show();
	}

	public ImagePlus getControlImage() {
		return getControlImage(image, new Point4f((float)x0, (float)y0, (float)z0, (float)r));
	}

	public static ImagePlus getControlImage(ImagePlus image, Point4f model) {
		Roi roi = image.getRoi();
		image.killRoi();
		ImagePlus imp = new Duplicator().run(image);
		ImageConverter.setDoScaling(true);
		new StackConverter(imp).convertToGray8();
		// imp = NaiveResampler.resample(imp, 4);

		IndexedTriangleMesh mesh = new IndexedTriangleMesh(
				MeshMaker.createSphere(model.x, model.y, model.z, model.w, 48, 48));
		imp = mesh.createOverlay(imp, 0xff0000);
		image.setRoi(roi);
		return imp;
	}

	public static ImagePlus getControlImage(ImagePlus image, Point4f model1, Point4f model2) {
		Roi roi = image.getRoi();
		image.killRoi();
		ImagePlus imp = new Duplicator().run(image);
		ImageConverter.setDoScaling(true);
		new StackConverter(imp).convertToGray8();
		// imp = NaiveResampler.resample(imp, 4);

		IndexedTriangleMesh mesh = new IndexedTriangleMesh(
				MeshMaker.createSphere(model1.x, model1.y, model1.z, model1.w, 48, 48));
		imp = mesh.createOverlay(imp, 0xff0000);

		mesh = new IndexedTriangleMesh(
				MeshMaker.createSphere(model2.x, model2.y, model2.z, model2.w, 48, 48));
		imp = mesh.createOverlay(imp, 0x00ff00);
		image.setRoi(roi);
		return imp;
	}


	public void fitMixture(double threshold) {
		int w = image.getWidth();
		int h = image.getHeight();
		int d = image.getStackSize();

		Calibration cal = image.getCalibration();
		double pw = cal.pixelWidth;
		double ph = cal.pixelHeight;
		double pd = cal.pixelDepth;

		ArrayList<Point3f> candidates = new ArrayList<Point3f>();
		double[] cog = new double[3];

		for(int zi = 0; zi < d; zi++) {
			ImageProcessor ip = image.getStack().getProcessor(zi + 1);
			float z = (float)(zi * pd);
			for(int yi = 0; yi < h; yi++) {
				float y = (float)(yi * ph);
				for(int xi = 0; xi < w; xi++) {
					if(ip.getf(xi, yi) < threshold)
						continue;
					if(image.getRoi() != null && !image.getRoi().contains(xi, yi))
						continue;
					float x = (float)(xi * pw);
					candidates.add(new Point3f(x, y, z));
					cog[0] += x;
					cog[1] += y;
					cog[2] += z;
				}
			}
		}
		cog[0] /= candidates.size();
		cog[1] /= candidates.size();
		cog[2] /= candidates.size();

		float cx = w * (float)pw / 2;
		float cy = h * (float)ph / 2;
		float cz = d * (float)pd / 2;

		cx = (float)cog[0];
		cy = (float)cog[1];
		cz = (float)cog[2];

		Point4f model1 = new Point4f(cx, cy, cz, 350);
		Point4f model2 = new Point4f(cx, cy, cz, 550);
		fitMixture(candidates, model1, model2);
		// getControlImage(image, model1, model2).show();

		Point4f smaller = model2.w < model1.w ? model2 : model1;

		x0 = smaller.x;
		y0 = smaller.y;
		z0 = smaller.z;
		r = smaller.w;
	}

	public void fitRansac(double threshold) {
		int w = image.getWidth();
		int h = image.getHeight();
		int d = image.getStackSize();

		Calibration cal = image.getCalibration();
		double pw = cal.pixelWidth;
		double ph = cal.pixelHeight;
		double pd = cal.pixelDepth;

		ArrayList<Point3f> candidates = new ArrayList<Point3f>();

		for(int zi = 0; zi < d; zi++) {
			ImageProcessor ip = image.getStack().getProcessor(zi + 1);
			float z = (float)(zi * pd);
			for(int yi = 0; yi < h; yi++) {
				float y = (float)(yi * ph);
				for(int xi = 0; xi < w; xi++) {
					if(ip.getf(xi, yi) < threshold)
						continue;
					if(image.getRoi() != null && !image.getRoi().contains(xi, yi))
						continue;
					float x = (float)(xi * pw);
					candidates.add(new Point3f(x, y, z));
				}
			}
		}

		int minNumMatches = candidates.size() / 50;
		int minNumInliers = candidates.size() / 20;
		int iterations = 20;
		double epsilon = 10;
		ransac(candidates,
				new ArrayList<Point3f>(),
				minNumMatches,
				iterations,
				epsilon,
				minNumInliers);

	}

	public void fit(double threshold) {
		int w = image.getWidth();
		int h = image.getHeight();
		int d = image.getStackSize();

		Calibration cal = image.getCalibration();
		double pw = cal.pixelWidth;
		double ph = cal.pixelHeight;
		double pd = cal.pixelDepth;

		double[][] XX = new double[4][4];
		double[][] XY = new double[4][1];

		for(int zi = 0; zi < d; zi++) {
			ImageProcessor ip = image.getStack().getProcessor(zi + 1);
			double z = zi * pd;
			for(int yi = 0; yi < h; yi++) {
				double y = yi * ph;
				for(int xi = 0; xi < w; xi++) {
					if(ip.getf(xi, yi) < threshold)
						continue;
					if(image.getRoi() != null && !image.getRoi().contains(xi, yi))
						continue;
					double x = xi * pw;

					double xx = x * x;
					double yy = y * y;
					double zz = z * z;
					double xyz = -xx - yy - zz;
					XX[0][0] += xx;
					XX[0][1] += x * y;
					XX[0][2] += x * z;
					XX[0][3] += x;
					XX[1][1] += yy;
					XX[1][2] += y * z;
					XX[1][3] += y;
					XX[2][2] += z * z;
					XX[2][3] += z;
					XX[3][3] += 1;
					XY[0][0] += x * xyz;
					XY[1][0] += y * xyz;
					XY[2][0] += z * xyz;
					XY[3][0] += xyz;
				}
			}
		}
		XX[1][0] = XX[0][1];
		XX[2][0] = XX[0][2];
		XX[3][0] = XX[0][3];
		XX[2][1] = XX[1][2];
		XX[3][1] = XX[1][3];
		XX[3][2] = XX[2][3];

		int n = (int)XX[3][3];
		for(int r = 0; r < 4; r++)
			for(int c = 0; c < 4; c++)
				XX[r][c] /= n;
		for(int r = 0; r < 4; r++)
			XY[r][0] /= n;

		Matrix m = new Matrix(XX);
		double[][] s = m.inverse().times(new Matrix(XY)).getArray();

		double a = s[0][0], b = s[1][0], c = s[2][0], e = s[3][0];

		x0 = -a / 2;
		y0 = -b / 2;
		z0 = -c / 2;
		r = Math.sqrt(((a * a + b * b + c * c) / 4) - e);
	}

	public final boolean ransac(
			ArrayList<Point3f> candidates,
			Collection<Point3f> inliers,
			int minNumMatches,
			int iterations,
			double epsilon,
			int minNumInliers ) {

		if(candidates.size() < minNumMatches)
			throw new RuntimeException("At least " +
					minNumMatches + " data points required." );

		inliers.clear();

		float bestMSE = Float.MAX_VALUE;

		for(int i = 0; i < iterations; i++) {
			ArrayList<Point3f> minMatches = new ArrayList<Point3f>();
			HashSet<Point3f> remaining = new HashSet<Point3f>();
			remaining.addAll(candidates);

			// choose minNumMatches disjunctive matches randomly
			for(int j = 0; j < minNumMatches; j++) {
				Point3f p = null;
				do {
					p = candidates.get((int)(Math.random()
							* candidates.size()));
				} while(!remaining.contains(p));
				minMatches.add(p);
				remaining.remove(p);
			}

			Point4f model = new Point4f();

			// fit this first matches
			fit(minMatches, model);
			System.out.println(model);

			// add all remaining candidate matches which are
			// explained by the model to the consensus set
			for(Point3f p : remaining) {
				if(distance(p, model) <= epsilon)
					minMatches.add(p);
			}

			System.out.println("minMatches.size() = " + minMatches.size());
			System.out.println("minNumInliers = " + minNumInliers);
			// continue if we found a bad model
			if(minMatches.size() < minNumInliers)
				continue;

			// we found a potentially good model
			fit(minMatches, model);
			System.out.println("refined: " + model);
			float mse = (float)mse(minMatches, model);
			if(mse < bestMSE) {
				inliers.clear();
				inliers.addAll(minMatches);
				bestMSE = mse;
				this.x0 = model.x;
				this.y0 = model.y;
				this.z0 = model.z;
				this.r = model.w;
			}
		}
		return !inliers.isEmpty();
	}

	public static double distance(Point3f p, Point4f model) {
		double dx = p.x - model.x;
		double dy = p.y - model.y;
		double dz = p.z - model.z;
		return Math.abs(Math.sqrt(dx * dx + dy * dy + dz * dz) - model.w);
	}

	public static double mse(Iterable<Point3f> pts, Point4f model) {
		double mse = 0;
		int size = 0;
		for(Point3f p : pts) {
			double d = distance(p, model);
			mse += d * d;
			size++;
		}
		return mse / size;
	}

	public static void fit(Iterable<Point3f> pts, Point4f model) {
		double[][] XX = new double[4][4];
		double[][] XY = new double[4][1];

		for(Point3f p : pts) {
			double z = p.z;
			double y = p.y;
			double x = p.x;

			double xx = x * x;
			double yy = y * y;
			double zz = z * z;
			double xyz = -xx - yy - zz;
			XX[0][0] += xx;
			XX[0][1] += x * y;
			XX[0][2] += x * z;
			XX[0][3] += x;
			XX[1][1] += yy;
			XX[1][2] += y * z;
			XX[1][3] += y;
			XX[2][2] += z * z;
			XX[2][3] += z;
			XX[3][3] += 1;
			XY[0][0] += x * xyz;
			XY[1][0] += y * xyz;
			XY[2][0] += z * xyz;
			XY[3][0] += xyz;
		}
		XX[1][0] = XX[0][1];
		XX[2][0] = XX[0][2];
		XX[3][0] = XX[0][3];
		XX[2][1] = XX[1][2];
		XX[3][1] = XX[1][3];
		XX[3][2] = XX[2][3];

		int n = (int)XX[3][3];
		for(int r = 0; r < 4; r++)
			for(int c = 0; c < 4; c++)
				XX[r][c] /= n;
		for(int r = 0; r < 4; r++)
			XY[r][0] /= n;

		Matrix m = new Matrix(XX);
		double[][] s = m.inverse().times(new Matrix(XY)).getArray();

		double a = s[0][0], b = s[1][0], c = s[2][0], e = s[3][0];

		model.x = (float)(-a / 2);
		model.y = (float)(-b / 2);
		model.z = (float)(-c / 2);
		model.w = (float)Math.sqrt(((a * a + b * b + c * c) / 4) - e);
	}

	public void fitMixture(Iterable<Point3f> pts, Point4f model1, Point4f model2) {
		ArrayList<Point3f> pts1 = new ArrayList<Point3f>();
		ArrayList<Point3f> pts2 = new ArrayList<Point3f>();
		double prevmse = Float.MAX_VALUE;
		while(true) {
			pts1.clear();
			pts2.clear();

			// calculate assignments
			for(Point3f p : pts) {
				double d1 = distance(p, model1);
				double d2 = distance(p, model2);
				if(Math.min(d1, d2) > 50)
					continue;
				if(d1 < d2)
					pts1.add(p);
				else
					pts2.add(p);
			}
			// update models
			fit(pts1, model1);
			fit(pts2, model2);

			double mse = mse(pts1, model1) + mse(pts2, model2);
			if(mse >= prevmse)
				break;
			prevmse = mse;
		}
	}
}
