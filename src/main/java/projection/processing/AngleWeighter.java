package projection.processing;

import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import javax.vecmath.Point3f;

public class AngleWeighter implements FusionWeight {

	public static final int X_AXIS = 0;
	public static final int Y_AXIS = 1;
	public static final int Z_AXIS = 2;

	private final int aperture;
	private final Point3f center;
	private final double angle;
	private final int axis;
	public static final double overlap = 20;
	public static final double overlap2 = overlap / 2;

	public AngleWeighter(int axis, double angle, int aperture, Point3f center) {
		this.axis = axis;
		this.angle = angle;
		this.aperture = aperture;
		this.center = center;
	}

	public double getAngle(double dx, double dy, double dz) {
		double d1, d2;
		switch(axis) {
		case Z_AXIS: d1 =  dx; d2 =  dy; break;
		case Y_AXIS: d1 =  dx; d2 = -dz; break;
		case X_AXIS: d1 = -dz; d2 =  dy; break;
		default: throw new IllegalArgumentException();
		}
		double a = 180.0 * Math.atan2(d2, d1) / Math.PI;

		a -= angle;
		if(a > 180)
			return a - 360;
		if(a < -180)
			return a + 360;
		return a;
	}

	public boolean inLowerOverlap(double angle) {
		return angle > -aperture/2.0 - overlap2 && angle < -aperture/2.0 + overlap2;
	}

	public boolean inUpperOverlap(double angle) {
		return angle > aperture/2.0 - overlap2 && angle < aperture/2.0 + overlap2;
	}


	/*
	 *          Camera
	 *            |
	 *            |
	 *
	 *            x
	 *         _________
	 *      y /        /|
	 *       /        / |
	 *      /________/  |
	 *      |        |  /
	 *    z |        | /
	 *      |________|/
	 */
	@Override
	public float getWeight(float px, float py, float pz) {

		double dx = px - center.x;
		double dy = py - center.y;
		double dz = pz - center.z;

		double angle = getAngle(dx, dy, dz);

		return getWeight(angle);
	}

	public float getWeight(double angle) {
		// inside
		if(angle > -aperture/2.0 + overlap2 && angle < aperture/2.0 - overlap2)
			return 1f;

		// outside
		if(angle < -aperture/2.0 - overlap2 || angle > aperture/2.0 + overlap2)
			return 0f;


		// within the overlap
		double eps = 0.0001;
		double k = 1 / overlap2 * Math.log(1.0 / eps - 1.0);
		if(angle < 0) {
			double t = angle + aperture/2.0;
			return (float)(1.0 / (1 + Math.exp(-k * t)));
		}
		double t = angle - aperture/2.0;
		return 1f - (float)(1.0 / (1 + Math.exp(-k * t)));
	}

	public void weightImage(ImagePlus image) {
		double pd = image.getCalibration().pixelDepth;
		double pw = image.getCalibration().pixelWidth;
		double ph = image.getCalibration().pixelHeight;
		int w = image.getWidth(), h = image.getHeight(), d = image.getStackSize();
		for(int z = 0; z < d; z++) {
			float rz = (float)(z * pd);
			ImageProcessor ip = image.getStack().getProcessor(z + 1);
			for(int y = 0; y < h; y++) {
				float ry = (float)(y * ph);
				for(int x = 0; x < w; x++) {
					float rx = (float)(x * pw);
					float v = ip.getf(x, y);
					ip.setf(x,  y, v * getWeight(rx, ry, rz));
				}
			}
		}
	}

	public static ImagePlus sum(ImagePlus[] images) {
		int w = images[0].getWidth(), h = images[0].getHeight();
		int d = images[0].getStackSize();
		ImagePlus res = IJ.createImage("sum", w, h, d, images[0].getBitDepth());
		int wh = w * h;
		for(ImagePlus imp : images) {
			for(int z = 0; z < d; z++) {
				ImageProcessor ip = imp.getStack().getProcessor(z + 1);
				ImageProcessor out = res.getStack().getProcessor(z + 1);
				for(int i = 0; i < wh; i++)
					out.setf(i, out.getf(i) + ip.getf(i));
			}
		}
		res.setCalibration(images[0].getCalibration().copy());
		return res;
	}

	public static void main(String[] args) {
		new ImageJ();
		int[] angles = new int[] {-135, -45, 45, 135};
		for(int a = 0; a < angles.length; a++) {
			AngleWeighter aw = new AngleWeighter(X_AXIS, angles[a], 45, new Point3f(50, 50, 50));
			int w = 100, h = 100, d = 100;
			ImageStack stack = new ImageStack(w, h);
			for(int z = 0; z < d; z++) {
				ImageProcessor p = new FloatProcessor(w, h);
				for(int y = 0; y < h; y++) {
					for(int x = 0; x < w; x++) {
						p.setf(x, y, aw.getWeight(x, y, z));
					}
				}
				stack.addSlice("", p);
			}
			new ImagePlus("aw_" + angles[a], stack).show();
		}

		AngleWeighter aw = new AngleWeighter(X_AXIS, 135, 90, new Point3f(50, 50, 50));
		double[] xd = new double[360];
		double[] yd = new double[360];
		for(int i = 0; i < xd.length; i++) {
			xd[i] = i - 180;
			yd[i] = aw.getWeight(xd[i]);
		}
		Plot p = new Plot("weights", "angle", "weight", xd, yd);
		p.show();

		int nAngles = 2;
		int angleInc = 45;
		int aperture = 90 / nAngles;
		int CAMERA1 = 0, CAMERA2 = 1;
		int LEFT = 0, RIGHT = 1;
		Point3f center = new Point3f(50, 50, 50);
		AngleWeighter[][][] weights = new AngleWeighter[2][2][nAngles];
		for(int a = 0; a < nAngles; a++) {
			weights[CAMERA1][LEFT] [a] = new AngleWeighter(AngleWeighter.X_AXIS,  135 + a * angleInc, aperture, center);
			weights[CAMERA1][RIGHT][a] = new AngleWeighter(AngleWeighter.X_AXIS, -135 + a * angleInc, aperture, center);
			weights[CAMERA2][LEFT] [a] = new AngleWeighter(AngleWeighter.X_AXIS,   45 + a * angleInc, aperture, center);
			weights[CAMERA2][RIGHT][a] = new AngleWeighter(AngleWeighter.X_AXIS,  -45 + a * angleInc, aperture, center);
		}
		int w = 1, h = 100, d = 100;
		for(int a = 0; a < 2; a++) {
			for(int cam = 0; cam < 2; cam++) {
				for(int ill = 0; ill < 2; ill++) {
					ImageStack stack = new ImageStack(w, h);
					for(int z = 0; z < d; z++) {
						ImageProcessor ip = new FloatProcessor(w, h);
						for(int y = 0; y < h; y++) {
							for(int x = 0; x < w; x++) {
								ip.setf(x, y, weights[cam][ill][a].getWeight(x, y, z));
							}
						}
						stack.addSlice("", ip);
					}
					new ImagePlus("angle" + a + "_ill" + ill + "_cam" + cam, stack).show();
				}
			}
		}
	}
}
