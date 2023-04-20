package projection.processing;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.gui.WaitForUserDialog;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import ij3d.Image3DUniverse;

import java.awt.geom.GeneralPath;
import java.awt.geom.PathIterator;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;

import org.scijava.java3d.Transform3D;
import org.scijava.vecmath.Matrix4f;
import org.scijava.vecmath.Point2f;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3f;

import projection.util.TransformIO;
import projection.viz.CustomContent;
import projection.viz.SphereProjectionViewer;

import com.jhlabs.map.proj.AugustProjection;
import com.jhlabs.map.proj.BonneProjection;
import com.jhlabs.map.proj.GallProjection;
import com.jhlabs.map.proj.Kavraisky7Projection;
import com.jhlabs.map.proj.MercatorProjection;
import com.jhlabs.map.proj.OrthographicAzimuthalProjection;
import com.jhlabs.map.proj.WinkelTripelProjection;


public class Map_Projection extends TimelapseProcessor implements PlugIn {

	public static final int MERCATOR              = 0;
	public static final int GALLPETER             = 1;
	public static final int KAVRAYSKIY            = 2;
	public static final int WINKEL_TRIPEL         = 3;
	public static final int AUGUSTUS_EPICYCLOIDAL = 4;
	public static final int BONNE                 = 5;
	public static final int ORTHO_AZIMUTAL        = 6;
	public static final int FULLER                = 7;

	public static final String[] MAP_TYPES = new String[] {
		"Mercator",
		"Gall-Peter",
		"Kavrayskiy",
		"Winkel Tripel",
		"Augustus Epicycloidal",
		"Bonne",
		"Orthogonal Azimutal",
		"Fuller"};

	public static String DEFAULT_DATA_DIR = "";

	@Override
	public void run(String arg) {
		GenericDialogPlus gd = new GenericDialogPlus("Create 2D Maps");
		gd.addDirectoryField("Data directory", DEFAULT_DATA_DIR);
		gd.addNumericField("#Layers", 1, 0);
		gd.addCheckbox("Create_coastlines", false);
		gd.addCheckbox("Create_Longitude/Latitude lines", true);
		gd.addCheckbox("Create_camera contribution overlay", true);
		gd.addCheckbox("Interactively choose initial transform", false);
		gd.addNumericField("Target width", 800, 0);
		for(int i = 0; i < MAP_TYPES.length; i++)
			gd.addCheckbox(MAP_TYPES[i], false);
		gd.showDialog();
		if(gd.wasCanceled())
			return;
		File datadir = new File(gd.getNextString());
		int nLayers = (int)gd.getNextNumber();
		boolean doCoast = gd.getNextBoolean();
		boolean doLines = gd.getNextBoolean();
		boolean doContributions = gd.getNextBoolean();
		boolean chooseInitial = gd.getNextBoolean();
		int tgtWidth = (int)gd.getNextNumber();

		if(!datadir.isDirectory()) {
			IJ.error(datadir + " is not a directory");
			return;
		}

		File objfile = new File(datadir, "Sphere.obj");
		if(!objfile.exists()) {
			IJ.error("Cannot find " + objfile.getAbsolutePath());
			return;
		}

		Matrix4f initial = new Matrix4f();
		initial.setIdentity();

		if(chooseInitial) {
			File previewdir = new File(datadir, "resampled");
			if(!previewdir.exists() || !previewdir.isDirectory())
				previewdir = datadir;

			Image3DUniverse univ = SphereProjectionViewer.show(previewdir.getAbsolutePath() + "/Sphere.obj", previewdir.getAbsolutePath(), null);
			new WaitForUserDialog("",
				"Please rotate the sphere to the desired orientation, then click OK").show();
			Transform3D trans = new Transform3D();
			CustomContent cc = (CustomContent)univ.getContent("bla");
			cc.getLocalRotate(trans);
			trans.get(initial);
			univ.close();
		} else if(new File(datadir, "initial_map_transform.mat").exists()) {
			try {
				initial = TransformIO.loadTransformation(new File(datadir, "initial_map_transform.mat"));
			} catch(Exception e) {
				e.printStackTrace();
			}
		}

		try {
			TransformIO.saveTransformation(new File(datadir, "initial_map_transform.mat"), initial);
		} catch(Exception e) {
			e.printStackTrace();
		}

		for(int i = 0; i < MAP_TYPES.length; i++) {
			if(!gd.getNextBoolean())
				continue;
			File outputdir = new File(datadir, MAP_TYPES[i]);
			if(outputdir.isDirectory()) {
				boolean cancelled = !IJ.showMessageWithCancel("Overwrite",
						outputdir + " already exists. Overwrite?");
				if(cancelled)
					return;
			} else {
				outputdir.mkdir();
			}

			try {
				prepareForProjection(nLayers, datadir, initial, i, outputdir, doCoast, doLines, doContributions, tgtWidth);
				createProjections(nLayers);
			} catch(Exception e) {
				IJ.error(e.getMessage());
				e.printStackTrace();
			}
		}
	}

	private SphericalMaxProjection smp;
	private GeneralProjProjection proj;
	private File datadir, outputdir;
	private int nVertices;
	private boolean doContributions, doLines, doCoast;
	private Matrix4f initial;
	private int nLayers;
	private FindMaxima findMaxima;

	public void prepareForProjection(
			final int nLayers,
			final File datadir,
			final Matrix4f initial,
			final int maptype,
			final File outputdir,
			final boolean doCoast,
			final boolean doLines,
			final boolean doContributions,
			final int tgtWidth) throws IOException {

		this.smp = new SphericalMaxProjection(new File(datadir, "Sphere.obj").getAbsolutePath(), initial);
		this.nLayers = nLayers;
		this.nVertices = smp.getSphere().nVertices;
		this.doLines = doLines;
		this.doCoast = doCoast;
		this.doContributions = doContributions;
		this.initial = initial;
		this.datadir = datadir;
		this.outputdir = outputdir;
		if(!outputdir.exists())
			outputdir.mkdirs();
		if(doContributions)
			new File(outputdir, "contributions").mkdir();
		if(doLines)
			new File(outputdir, "lines").mkdir();
		if(doCoast)
			new File(outputdir, "coast").mkdir();

		switch(maptype) {
			case MERCATOR:              proj = new GeneralProjProjection(new MercatorProjection());   break;
			case GALLPETER:             proj = new GeneralProjProjection(new GallProjection()); break;
			case KAVRAYSKIY:            proj = new GeneralProjProjection(new Kavraisky7Projection()); break;
			case WINKEL_TRIPEL:         proj = new GeneralProjProjection(new WinkelTripelProjection()); break;
			case AUGUSTUS_EPICYCLOIDAL: proj = new GeneralProjProjection(new AugustProjection()); break;
			case BONNE:                 proj = new GeneralProjProjection(new BonneProjection()); break;
			case ORTHO_AZIMUTAL:        proj = new GeneralProjProjection(new OrthographicAzimuthalProjection()); break;
			case FULLER:                proj = new FullerProjection(); break;
			default: throw new IllegalArgumentException("Unsupported map type: " + maptype);
		}
		proj.prepareForProjection(smp, tgtWidth);
		findMaxima = new FindMaxima(smp);
	}

	ArrayList<Point3f> readPoints(File f) throws IOException {
		ArrayList<Point3f> list = new ArrayList<Point3f>();
		BufferedReader in = new BufferedReader(new FileReader(f));
		String line = null;
		while((line = in.readLine()) != null) {
			String[] toks = line.split(",");
			list.add(new Point3f(
					(float)Double.parseDouble(toks[0].trim()),
					(float)Double.parseDouble(toks[0].trim()),
					(float)Double.parseDouble(toks[0].trim())));
		}
		in.close();
		return list;
	}

	@Override
	public void processTimepoint(int tp) throws IOException {
		System.out.println("fuller: time point " + tp);
		String basename = String.format("tp%04d", tp);

		File matdir = new File(datadir, "transformations");
		File contribindir = new File(datadir, "contributions");
		File contriboutdir = new File(outputdir, "contributions");
		File linesdir = new File(outputdir, "lines");
		File coastdir = new File(outputdir, "coast");

		File matfile = new File(matdir, basename + ".matrix");
		File contribin = new File(contribindir, basename + ".vertices");
		if(!contribin.exists())
			contribin = new File(datadir, "contributions.vertices");

		int w = 0, h = 0;
		for(int l = 0; l < nLayers; l++) {
			File outfile = new File(outputdir, String.format("tp%04d_%02d.tif", tp, l));
			File infile = new File(datadir, String.format("tp%04d_%02d.vertices", tp, l));

			short[] maxima = SphericalMaxProjection.loadShortData(infile.getAbsolutePath(), nVertices);

			Overlay overlay = null;
			File ptsfile = new File(datadir, String.format("tp%04d_%02d.pts", tp, l));
			if(ptsfile.exists()) {
				ArrayList<Point3f> pts = readPoints(ptsfile);
				short[][] vals = new short[1][smp.getSphere().nVertices];
				vals[0] = maxima;
				// ArrayList<Point3f> pts = Register_.getPoints(findMaxima, smp, vals);
				float[] xs = new float[pts.size()];
				float[] ys = new float[pts.size()];
				Point2f polar = new Point2f();
				double[] p2d = new double[2];
				int i = 0;
				for(Point3f p : pts) {
					smp.getPolar(p, polar);
					toDeg(polar);
					p2d[0] = polar.x;
					p2d[1] = polar.y;
					proj.transform(p2d);
					xs[i] = (float)p2d[0];
					ys[i] = (float)p2d[1];
					i++;
				}
				PointRoi roi = new PointRoi(xs, ys);
				overlay = new Overlay(roi);
			}

			// maxima = smp.applyTransform(get90DegRot(smp), maxima);
			ImageProcessor ip = proj.project(maxima);
			w = ip.getWidth();
			h = ip.getHeight();
			ImagePlus image = new ImagePlus("", ip);
			image.setOverlay(overlay);

			// coast = transform(smp, mat, coast);
			IJ.save(image, outfile.getAbsolutePath());
		}
		if(doLines) {
			GeneralPath lines = GeneralProjProjection.createLines();
			Matrix4f mat = new Matrix4f(initial);
			if(matfile.exists())
				mat.mul(TransformIO.loadTransformation(matfile));
			else
				System.out.println(matfile + " does not exist");
			lines = transform(smp, mat, lines);

			// rotate the lines by 90 deg
			Matrix4f rot = get90DegRot(smp);
			lines = transform(smp, rot, lines);

			lines = proj.transform(lines);
			GeneralProjProjection.savePath(lines, new File(linesdir, basename + ".eps").getAbsolutePath(), proj.getWidth(), proj.getHeight(), false);
			ImageProcessor lip = new ByteProcessor(w, h);
			proj.drawInto(lip, 255, 1, lines);
			IJ.save(new ImagePlus("", lip), new File(linesdir, basename + ".tif").getAbsolutePath());
		}
		if(doCoast) {
			URL url = getClass().getResource("coast.dat");
			GeneralPath coast = GeneralProjProjection.readDatFile(url);
			Matrix4f mat = new Matrix4f(initial);
			if(matfile.exists())
				mat.mul(TransformIO.loadTransformation(matfile));
			else
				System.out.println(matfile + " does not exist");
			coast = transform(smp, mat, coast);

			// rotate the lines by 90 deg
			Matrix4f rot = get90DegRot(smp);
			coast = transform(smp, rot, coast);

			coast = proj.transform(coast);
			GeneralProjProjection.savePath(coast, new File(coastdir, basename + ".eps").getAbsolutePath(), proj.getWidth(), proj.getHeight(), false);
		}
		if(doContributions) {
			int[] contribs = SphericalMaxProjection.loadIntData(contribin.getAbsolutePath(), nVertices);
			ImageProcessor overlay = proj.projectColor(contribs);
			IJ.save(new ImagePlus("", overlay), new File(contriboutdir, basename + ".tif").getAbsolutePath());
		}
	}

	public void createProjections(final int nLayers) {
		final String[] files = datadir.list(new FilenameFilter() {
			@Override
			public boolean accept(File dir, String name) {
				return name.startsWith("tp") && name.endsWith("_00.vertices");
			}

		});
		processParallel(files.length);
	}

	public static Matrix4f get90DegRot(SphericalMaxProjection smp) {
		Matrix4f rot = new Matrix4f();
		Matrix4f cen = new Matrix4f();

		rot.rotZ((float)(Math.PI / 2));
		cen.set(new Vector3f(smp.center));

		rot.mul(cen, rot);
		cen.invert();
		rot.mul(rot, cen);

		return rot;
	}

	/**
	 * In- and output coordinates are in degrees.
	 * @param smp
	 * @param trans
	 * @param path
	 * @return
	 */
	public static GeneralPath transform(SphericalMaxProjection smp, Matrix4f trans, GeneralPath path) {
		GeneralPath out = new GeneralPath();
		PathIterator it = path.getPathIterator(null);
		float[] seg = new float[6];

		Point2f pin = new Point2f();

		while(!it.isDone()) {
			int l = it.currentSegment(seg);
			pin.x = seg[0];
			pin.y = seg[1];
			toRad(pin);
			transform(smp, trans, pin);
			toDeg(pin);

			if(l == PathIterator.SEG_MOVETO)
				out.moveTo(pin.x, pin.y);
			else
				out.lineTo(pin.x, pin.y);
			it.next();
		}
		return out;
	}

	public static void toRad(Point2f polar) {
		polar.x = (float)(polar.x / 180 * Math.PI);
		polar.y = (float)(polar.y / 180 * Math.PI);
	}

	public static void toDeg(Point2f polar) {
		polar.x = (float)(polar.x / Math.PI * 180);
		polar.y = (float)(polar.y / Math.PI * 180);
	}

	public static void transform(SphericalMaxProjection smp, Matrix4f trans, Point2f polar) {
		Point3f coord = new Point3f();
		Point3f center = smp.getCenter();
		float radius = smp.getRadius();
		getPoint(center, radius, polar.x, polar.y, coord);
		trans.transform(coord);
		smp.getPolar(coord, polar);
	}

	public static void getPoint(Point3f center, float radius, float longitude, float latitude, Point3f ret) {
		double sinLong = Math.sin(longitude),
			cosLong = Math.cos(longitude),
			sinLat = Math.sin(latitude),
			cosLat = Math.cos(latitude);
		ret.z = (float)(center.z + radius * cosLat * cosLong);
		ret.x = (float)(center.x - radius * cosLat * sinLong);
		ret.y = (float)(center.y - radius * sinLat);
	}
}
