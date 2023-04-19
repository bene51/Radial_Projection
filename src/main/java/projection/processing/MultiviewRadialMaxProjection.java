package projection.processing;


import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import org.scijava.vecmath.Matrix4f;
import org.scijava.vecmath.Point3f;

import projection.util.TransformIO;


public class MultiviewRadialMaxProjection {

	protected final String outputdir;
	protected final Configuration conf;
	// protected final int nChannels;
	protected final int[] channels;

	protected SphericalMaxProjection[] smp;
	protected MaximumProjector[] projectors;

	private static final int[] createIncrementingArray(int n) {
		int[] ret = new int[n];
		for(int i = 0; i < n; i++) {
			ret[i] = i;
		}
		return ret;
	}

	public MultiviewRadialMaxProjection(String outputdir,
			Configuration conf,
			int nChannels) {
		this(outputdir, conf, createIncrementingArray(nChannels));
	}

	public MultiviewRadialMaxProjection(String outputdir,
			Configuration conf,
			int[] channels) {
		if(!outputdir.endsWith(File.separator))
			outputdir += File.separator;

		this.conf = conf;
		this.outputdir = outputdir;
		this.channels = channels;

		try {
			// initialize the maximum projections
			initSphericalMaximumProjection();
			// this.maxima = new short[conf.nAngles][conf.nLayers][smp[0].getSphere().nVertices];
		} catch(Exception e) {
			throw new RuntimeException("Cannot load transformations.", e);
		}

		// save the sphere geometry
		String spherepath = new File(outputdir, "Sphere.obj").getAbsolutePath();
		try {
			smp[0].saveSphere(spherepath);
			saveTransformations(conf.transformations);
			saveVertexIndices();
		} catch(Exception e) {
			throw new RuntimeException("Cannot save sphere: " + spherepath, e);
		}
	}

	private void saveVertexIndices() throws IOException {
		for(int c : channels) {
			File dir = new File(outputdir, "channel" + c);
			for(int a = 0; a < conf.angles.length; a++)
				projectors[a].saveVertexIndices(new File(dir, getIndicesFileName(conf.angleNames, a)));
		}
	}

	private void saveTransformations(Matrix4f[] transformations) throws IOException {
		for(int c : channels) {
			File dir = new File(outputdir, "channel" + c);
			System.out.println("make dir " + dir.getAbsolutePath());
			if(!dir.exists())
				dir.mkdir();
			for(int aIndex = 0; aIndex < conf.angles.length; aIndex++) {
				PrintStream out = new PrintStream(new File(dir, conf.angleNames[aIndex] + ".mat"));
				out.println(conf.angleNames[aIndex]);
				out.println(conf.angles[aIndex]);
				out.println(conf.apertures[aIndex]);
				out.println(TransformIO.toString(transformations[aIndex]));
				out.close();
			}
		}
	}

	/**
	 * Returns the output directory.
	 * @return
	 */
	public String getOutputDirectory() {
		return outputdir;
	}

	/**
	 * Process the next image, given in form of a short[] array
	 * @param ip image as short[] array.
	 */
	public void process(short[] ip, int tp, int aIndex, int channel, int z) {

		// Start of stack
		if(z == 0)
			projectors[aIndex].resetMaxima();

		// do the projection
		projectors[aIndex].projectPlaneMultilayer(z, ip);

		// Not end of stack: nothing else to do
		if(z < conf.d - 1)
			return;

		// save the result
		for(int l = 0; l < conf.nLayers; l++) {
			String filename = getFileName(tp, conf.angleNames, aIndex, l);
			File vpath = new File(outputdir, "channel" + channel);
			vpath = new File(vpath, filename);
			try {
				projectors[aIndex].saveVertices(l, vpath);
			} catch(Exception e) {
				throw new RuntimeException("Cannot save " + vpath);
			}
		}
	}

	public static String getFileName(int tp, String[] angleNames, int aIndex, int layer) {
		return String.format("tp%04d_%s_l%02d.vertices", tp, angleNames[aIndex], layer);
	}

	public static String getIndicesFileName(String[] angleNames, int aIndex) {
		return String.format("%s.indices", angleNames[aIndex]);
	}

	/**
	 * @param transform Array with one transformation for each angle;
	 * the first angle is not transformed, all following angles are transformed;
	 * example:
	 * 	2 samples,
	 * 	2 angles each
	 *
	 * called twice with 2 transforms:
	 * 	- [I], [a1->a0]
	 */
	private SphericalMaxProjection[] initSphericalMaximumProjection() throws IOException {
		int subd = (int)Math.round(conf.radius / (Math.min(conf.pw, Math.min(conf.ph, conf.pd))));
		int levels = (int)Math.floor(Math.log(subd) / Math.log(2));

		System.out.print("Create sphere");
		long start = System.currentTimeMillis();
		Point3f center = new Point3f((float)conf.cx, (float)conf.cy, (float)conf.cz);
		IndexedTriangleMesh sphere = createSphere(center, (float)conf.radius, levels);
		long end = System.currentTimeMillis();
		System.out.println("... took " + (end - start) + " ms: " + sphere.nVertices + " vtcs");

		smp = new SphericalMaxProjection[conf.nAngles];
		projectors = new MaximumProjector[conf.nAngles];

		for(int a = 0; a < conf.nAngles; a++) {
			Matrix4f transform = conf.transformations[a];
			Point3f cen = new Point3f(center);
			if(transform != null)
				transform.transform(cen);

			System.out.println("Create SMP");
			smp[a] = new SphericalMaxProjection(sphere, center, (float)conf.radius, transform, false);
			projectors[a] = new MaximumProjector(smp[a]);
			System.out.println("prepare for projection");
			projectors[a].prepareForProjection(
				conf.w, conf.h, conf.d,
				conf.pw, conf.ph, conf.pd,
				conf.layerwidth, conf.nLayers, conf.ratioInside,
				new AngleWeighter(AngleWeighter.X_AXIS, conf.angles[a], conf.apertures[a], cen));
		}
		return smp;
	}

	/**
	 * Create a spherical mesh using a buckyball approximation
	 * @param center The center of the sphere.
	 * @param radius The radius of the sphere.
	 * @param subd The number of subdivisions of the icosahedron.
	 * @return
	 */
	private static IndexedTriangleMesh createSphere(Point3f center, float radius, int levels) {
		// calculate the sphere coordinates
		Icosahedron icosa = new Icosahedron(radius);

		IndexedTriangleMesh sphere = icosa.createBuckyballFast(radius, levels);
		for(Point3f p : sphere.getVertices())
			p.add(center);
		return sphere;
	}
}

