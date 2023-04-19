package projection.processing;

import org.scijava.vecmath.Matrix4f;
import org.scijava.vecmath.Point3f;

import vib.FastMatrix;

public class HorizontalAligner {

	// initial: euler angles
	public static final Matrix4f align(Point3f[] vertices, Point3f center, int[] maxima, double[] initial, int threshold) {
		Matrix4f m = new Matrix4f();
		new HorizontalAligner(vertices, center, maxima, threshold).optimize(initial, m);
		return m;
	}

	public static final Matrix4f align(Point3f[] vertices, Point3f center, int[] maxima, int threshold) {
		return new HorizontalAligner(vertices, center, maxima, threshold).optimize();
	}

	private final Point3f[] vertices;
	private final int[] maxima;
	private final Point3f center;
	private final int threshold;

	public HorizontalAligner(Point3f[] vertices, Point3f center, int[] maxima, int threshold) {
		super();
		this.vertices = vertices;
		this.maxima = maxima;
		this.center = center;
		this.threshold = threshold;
	}

	/**
	 * Try all principal orientations
	 */
	public Matrix4f optimize() {
		double[][] initial = new double[24][];
		// Contruct initial matrices corresponding to the
		// 24 principal euler orientations
		float pi2 = (float)(Math.PI / 2.0);
		int n = 0;
		for(int i = 0; i < 4; i++)
			for(int j = 0; j < 4; j++)
				initial[n++] = new double[] { i * pi2, j * pi2, 0 };

		for(int i = 0; i < 4; i++)
			initial[n++] = new double[] { i * pi2, pi2, pi2 };

		for(int i = 0; i < 4; i++)
			initial[n++] = new double[] { i * pi2, 3 * pi2, pi2 };

		Matrix4f bestMatrix = new Matrix4f();
		double bestV = optimize(initial[0], bestMatrix);
		for(int i = 0; i < 24; i++) {
			Matrix4f m = new Matrix4f();
			double v = optimize(initial[i], m);
			if(v < bestV) {
				bestV = v;
				bestMatrix = m;
			}
		}
		return bestMatrix;
	}

	public double optimize(double[] initial, Matrix4f ret) {
		final double[] max = new double[] { Math.PI / 4, Math.PI / 4, Math.PI / 4 };

		double badness = new Optimizer() {
			@Override
			public double calculateDifference(double[] parameters) {
				return getDifference(getMatrix(parameters));
			}
		}.optimize(1.0, initial, max);

		ret.set(getMatrix(initial));
		return badness;
	}

	public Matrix4f getMatrix(double[] param) {
		FastMatrix fm = FastMatrix.rotateEulerAt(
			param[0], param[1], param[2],
			center.x, center.y, center.z);
		double[] d = fm.rowwise16();
		return new Matrix4f(
			(float)d[0], (float)d[1], (float)d[2], (float)d[3],
			(float)d[4], (float)d[5], (float)d[6], (float)d[7],
			(float)d[8], (float)d[9], (float)d[10], (float)d[11],
			(float)d[12], (float)d[13], (float)d[14], (float)d[15]);
	}

	public double getDifference(Matrix4f transform) {
		double sumY = 0;
		double sumY2 = 0;
		double total = 0;


		Point3f trans = new Point3f();
		for(int i = 0; i < vertices.length; i++) {
			Point3f v = vertices[i];
			transform.transform(v, trans);

			double m = maxima[i];
//			sumY  += trans.y * m;
//			sumY2 += trans.y * trans.y * m;
//			total += m;
			if(m > threshold) {
				sumY  += trans.y;
				sumY2 += trans.y * trans.y;
				total += 1;
			}
		}

		double Ey  = sumY / total;
		double Ey2 = sumY2 / total;

		double var = Ey2 - Ey * Ey;

		return var;
	}
}
