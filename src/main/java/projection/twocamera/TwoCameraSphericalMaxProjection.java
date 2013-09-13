package projection.twocamera;


import javax.vecmath.Matrix4f;
import javax.vecmath.Point3f;

import projection.processing.Configuration;
import projection.processing.MultiviewRadialMaxProjection;

/**
 * This class assumes we have two cameras and two angles
 * @author bschmid
 *
 */
public class TwoCameraSphericalMaxProjection extends MultiviewRadialMaxProjection {

	public static final int LEFT  = 0;
	public static final int RIGHT = 1;
	public static final int CAMERA1 = 0;
	public static final int CAMERA2 = 1;
	public static final int CAM1_ANGLE = 135;
	public static final int CAM2_ANGLE = 45;

	private static int[] getAngles(int nAngles, int camera) {
		int n = nAngles * 2; // 2 illumination sides
		int[] angles = new int[n];
		int i = 0;
		for(int a = 0; a < nAngles; a++) {
			angles[i++] = getAngle(camera, LEFT);
			angles[i++] = getAngle(camera, RIGHT);
		}
		return angles;
	}

	public static int getAngle(int camera, int ill) {
		int angle = camera == CAMERA1 ? CAM1_ANGLE : CAM2_ANGLE;
		if(ill == RIGHT)
			angle = -angle;
		return angle;
	}

	private static String[] getAngleNames(int nAngles, int angleInc, int camera) {
		int n = nAngles * 2; // 2 illumination sides
		String[] angles = new String[n];
		int i = 0;
		for(int a = 0; a < nAngles; a++) {
			angles[i++] = getAngleName(a, angleInc, camera, LEFT);
			angles[i++] = getAngleName(a, angleInc, camera, RIGHT);
		}
		return angles;
	}

	public static String getAngleName(int aIndex, int angleInc, int camera, int ill) {
		int angle =  getAngle(camera, ill);
		if(angle < 0)
			angle = 360 + angle;

		angle += aIndex * angleInc;
		return String.format("%04d", angle);
	}

	private static int[] getApertures(int nAngles) {
		int aperture = 90 / nAngles;
		int n = nAngles * 2; // 2 illumination sides
		int[] apertures = new int[n];
		for(int a = 0; a < n; a++)
			apertures[a] = aperture;
		return apertures;
	}

	private static Matrix4f[] makeMultiviewTransformations(Matrix4f[] transforms) {
		int n = transforms.length * 2; // 2 illumination sides
		Matrix4f[] ret = new Matrix4f[n];
		int i = 0;
		for(int a = 0; a < transforms.length; a++) {
			ret[i++] = transforms[a];
			ret[i++] = transforms[a];
		}
		return ret;
	}

	/**
	 *
	 * @param outputdir Output directory.
	 * @param camera Left or right camera.
	 * @param angleInc Angle increment.
	 * @param nAngles Number of angles.
	 * @param w Width of images taken.
	 * @param h Height of images taken.
	 * @param d Stack size.
	 * @param pw Pixel calibration.
	 * @param ph Pixel calibration.
	 * @param pd Pixel calibration.
	 * @param center Estimated center of the sphere.
	 * @param radius Estimated radius of the sphere
	 * @param transforms List of transformations, each specifying the forward transformation from angle[0] to angle[i].
	 */
	public TwoCameraSphericalMaxProjection(String outputdir,
			int camera,
			int angleInc, int nAngles,
			int w, int h, int d,
			double pw, double ph, double pd,
			Point3f center, float radius,
			double layerWidth,
			int nLayers,
			int nChannels,
			Matrix4f[] transforms) {

		super(outputdir, new Configuration(
			w, h, d,
			1, 2 * nAngles, nLayers, pw, ph, pd,
			center.x, center.y, center.z, radius, layerWidth,
			getAngles(nAngles, camera),
			getAngleNames(nAngles, angleInc, camera),
			makeMultiviewTransformations(transforms),
			getApertures(nAngles)),
			nChannels);
	}

	/**
	 * Process the next image, given in form of a short[] array
	 * @param ip image as short[] array.
	 */
	public void process(short[] ip, int tp, int aIndex, int channel, int z, int ill) {
		super.process(ip, tp, aIndex * 2 + ill, channel, z);
	}
}

