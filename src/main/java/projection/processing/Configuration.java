package projection.processing;


import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Properties;

import javax.vecmath.Matrix4f;

import projection.twocamera.TwoCameraSphericalMaxProjection;
import projection.util.TransformIO;

public class Configuration {

	enum names {
		width, height, nPlanes, nTimepoints, nAngles, nLayers, pixelwidth, pixelheight, pixeldepth, centerx, centery, centerz, radius, layerwidth,
	}

	public final int w, h, d, nTimepoints, nAngles, nLayers;
	public final double pw, ph, pd;
	public final double cx, cy, cz, radius;
	public final double layerwidth;

	public final int[] angles;
	public final String[] angleNames;
	public final Matrix4f[] transformations;
	public final int[] apertures;

	public Configuration(File f) throws Exception {
		Properties props = new Properties();
		FileReader fw = new FileReader(f);
		props.load(fw);
		this.w           = readInt(props, names.width.name());
		this.h           = readInt(props, names.height.name());
		this.d           = readInt(props, names.nPlanes.name());
		this.nTimepoints = readInt(props, names.nTimepoints.name());
		this.nAngles     = readInt(props, names.nAngles.name());
		this.nLayers     = readInt(props, names.nLayers.name());
		this.pw          = readDouble(props, names.pixelwidth.name());
		this.ph          = readDouble(props, names.pixelheight.name());
		this.pd          = readDouble(props, names.pixeldepth.name());
		this.cx          = readDouble(props, names.centerx.name());
		this.cy          = readDouble(props, names.centery.name());
		this.cz          = readDouble(props, names.centerz.name());
		this.radius      = readDouble(props, names.radius.name());
		this.layerwidth  = readDouble(props, names.layerwidth.name());

		angles = new int[nAngles];
		angleNames = new String[nAngles];
		transformations = new Matrix4f[nAngles];
		apertures = new int[nAngles];

		for(int a = 0; a < nAngles; a++) {
			angles[a] = readInt(props, "angle" + a + ".value");
			angleNames[a] = readString(props, "angle" + a + ".name");
			transformations[a] = TransformIO.fromString(readString(props, "angle" + a + ".transformation"));
			apertures[a] = readInt(props, "angle" + a + ".aperture");
		}
	}

	public Configuration(int w, int h, int d, int nTimepoints, int nAngles,
			int nLayers, double pw, double ph, double pd,
			double cx, double cy, double cz, double radius,
			double layerwidth, int[] angles, String[] angleNames,
			Matrix4f[] transformations, int[] apertures) {
		super();
		this.w = w;
		this.h = h;
		this.d = d;
		this.nTimepoints = nTimepoints;
		this.nAngles = nAngles;
		this.nLayers = nLayers;
		this.pw = pw;
		this.ph = ph;
		this.pd = pd;
		this.cx = cx;
		this.cy = cy;
		this.cz = cz;
		this.radius = radius;
		this.layerwidth = layerwidth;
		this.angles = angles;
		this.angleNames = angleNames;
		this.transformations = transformations;
		this.apertures = apertures;
	}

	private final String readString(Properties p, String key) throws Exception {
		if(!p.containsKey(key))
			throw new Exception("No '" + key + "' field in configuration file");
		return p.getProperty(key);
	}

	private final int readInt(Properties p, String key) throws Exception {
		String v = readString(p, key);
		try {
			return Integer.parseInt(v);
		} catch(Exception e) {
			throw new Exception("Cannot parse value for " + key + ": " + v);
		}
	}

	private final double readDouble(Properties p, String key) throws Exception {
		String v = readString(p, key);
		try {
			return Double.parseDouble(v);
		} catch(Exception e) {
			throw new Exception("Cannot parse value for " + key + ": " + v);
		}
	}

	public void saveConfiguration(File f) throws IOException {
		PrintStream out = new PrintStream(new FileOutputStream(f));
		out.println("#Created by RadialMaxProjection");
		out.println(names.width.name() + "=" + Integer.toString(w));
		out.println(names.height.name() + "=" + Integer.toString(h));
		out.println(names.nPlanes.name() + "=" + Integer.toString(d));
		out.println(names.nTimepoints.name() + "=" + Integer.toString(nTimepoints));
		out.println(names.nAngles.name() + "=" + Integer.toString(nAngles));
		out.println(names.nLayers.name() + "=" + Integer.toString(nLayers));
		out.println(names.pixelwidth.name() + "=" + Double.toString(pw));
		out.println(names.pixelheight.name() + "=" + Double.toString(ph));
		out.println(names.pixeldepth.name() + "=" + Double.toString(pd));
		out.println(names.centerx.name() + "=" + Double.toString(cx));
		out.println(names.centery.name() + "=" + Double.toString(cy));
		out.println(names.centerz.name() + "=" + Double.toString(cz));
		out.println(names.radius.name() + "=" + Double.toString(radius));
		out.println(names.layerwidth.name() + "=" + Double.toString(layerwidth));

		for(int a = 0; a < nAngles; a++)
			out.println("angle" + a + ".value" + "=" + Integer.toString(angles[a]));
		for(int a = 0; a < nAngles; a++)
			out.println("angle" + a + ".name" + "=" + angleNames[a]);
		for(int a = 0; a < nAngles; a++)
			out.println("angle" + a + ".transformation" + "=" + TransformIO.toString(transformations[a]));
		for(int a = 0; a < nAngles; a++)
			out.println("angle" + a + ".aperture" + "=" + Integer.toString(apertures[a]));

		out.close();
	}

	public static void main(String[] args) throws IOException {
		int w = 444, h = 440, d = 201;
		int nTimepoints = 1;
		int nAngles = 8;
		int nLayers = 1;
		double pw = 2.6, ph = 2.6, pd = 4;
		double cx = 587.82904, cy = 522.97534, cz = 367.6763, radius = 345.3111267089844;
		double layerwidth = 140;
		int[] apertures = new int[] {45, 45, 45, 45, 45, 45, 45, 45};
		int[] angles = new int[8];
		String[] angleNames = new String[8];
		Matrix4f[] transformations = new Matrix4f[8];

		for(int a = 0, i = 0; a < 2; a++) {
			Matrix4f trans = a == 0 ? null : new Matrix4f(
					0.9999492f, -0.004943217f, 0.008788793f, -4.8553543f,
					-0.0027186377f, 0.7071569f, 0.70705146f, -69.07996f,
					-0.009710164f, -0.7070394f, 0.7071075f, 499.9826f,
					0.0f, 0.0f, 0.0f, 1.0f);
			for(int cam = 0; cam < 2; cam++) {
				for(int ill = 0; ill < 2; ill++, i++) {
					angleNames[i] = TwoCameraSphericalMaxProjection.getAngleName(a, 45, cam, ill);
					System.out.println("angle " + a + " cam " + cam + " ill " + ill + ": " + angleNames[i]);
					angles[i] = TwoCameraSphericalMaxProjection.getAngle(cam, ill);
					transformations[i] = trans;
				}
			}
		}
		Configuration c = new Configuration(
				w, h, d, nTimepoints, nAngles, nLayers,
				pw, ph, pd, cx, cy, cz, radius, layerwidth,
				angles, angleNames, transformations, apertures);
		c.saveConfiguration(new File("/Users/bschmid/PostDoc/spim2_testdata/raw/RadialMaxProj.conf"));
	}
}
