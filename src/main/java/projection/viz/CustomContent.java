package projection.viz;

import ij.IJ;
import ij3d.Content;
import ij3d.ContentInstant;

import java.awt.Polygon;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.scijava.java3d.Canvas3D;
import org.scijava.java3d.GeometryArray;
import org.scijava.java3d.IndexedTriangleArray;
import org.scijava.java3d.Transform3D;
import org.scijava.java3d.TriangleArray;
import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point2d;
import org.scijava.vecmath.Point3d;
import org.scijava.vecmath.Point3f;

import projection.processing.SphericalMaxProjection;

import customnode.CustomMeshNode;
import customnode.CustomTriangleMesh;

public class CustomContent extends Content {

	private final CustomIndexedTriangleMesh mesh;
	public final File[] files;

	private float displayedMaximum = 0;
	private float displayedMinimum = 0;

	private float elevationFactor = 0.00001f;

	private final SphericalMaxProjection smp;
	private short[] maxima;

	private int currentIdx = 0;

	private boolean showMaxima = false;
	private boolean showColorOverlay = false;
	private byte[][] lut = null;

	private float maximaThreshold = 0;

	public CustomContent(String objfile, String vertexDir, String filenameContains) throws IOException {

		super("bla", 0);
		smp = new SphericalMaxProjection(objfile);

		List<File> tmp = new ArrayList<File>();
		tmp.addAll(Arrays.asList(new File(vertexDir).listFiles()));
		for(int i = tmp.size() - 1; i >= 0; i--) {
			String name = tmp.get(i).getName();
			if(!name.endsWith(".vertices"))
				tmp.remove(i);
			else if(filenameContains != null && !name.contains(filenameContains))
				tmp.remove(i);
		}
		this.files = new File[tmp.size()];
		tmp.toArray(files);
		Arrays.sort(files);

		int nVertices = smp.getSphere().nVertices;
		Color3f[] colors = new Color3f[nVertices];
		for(int i = 0; i < colors.length; i++)
			colors[i] = new Color3f(0, 1, 0);

		try {
			readColors(files[0]);
		} catch(Exception e) {
			throw new RuntimeException("Cannot load " + files[0], e);
		}
		mesh = new CustomIndexedTriangleMesh(smp.getSphere().getVertices(), colors, smp.getSphere().getFaces());
		CustomMeshNode node = new CustomMeshNode(mesh);

		ContentInstant content = getInstant(0);
		content.display(node);

		displayedMinimum = getCurrentMinimum();
		displayedMaximum = getCurrentMaximum();

		updateDisplayRange();
	}

	public SphericalMaxProjection getSMP() {
		return smp;
	}

	public Color3f getColor(int i) {
		return mesh.colors[i];
	}

	public void exportToPLY(String path) throws IOException {
		mesh.createNormalizedVersion(smp.getCenter(), smp.getRadius()).exportToPLY(path);
	}

	public void setElevationFactor(float f) {
		this.elevationFactor = f;
	}

	public float getElevationFactor() {
		return elevationFactor;
	}

	public void setLUT(byte[][] lut) {
		this.lut = lut;
		updateDisplayRange();
	}

	public String getCurrentFilePath() {
		return files[currentIdx].getAbsolutePath();
	}

	public File getCurrentFile() {
		return files[currentIdx];
	}

	public boolean areMaximaShown() {
		return showMaxima;
	}

	public void toggleShowMaxima() {
		showMaxima = !showMaxima;
		updateDisplayRange();
	}

	public boolean isShowAsColor() {
		return showColorOverlay;
	}

	public void toggleShowAsColor() {
		showColorOverlay = !showColorOverlay;
		updateDisplayRange();
	}

	public int[] getMaxima() {
		int[] m = new int[maxima.length];
		for(int i = 0; i < m.length; i++)
			m[i] = maxima[i] & 0xffff;
		return m;
	}

	public void smooth() {
		int nPasses = (int) IJ.getNumber("Passes", 5);
		smp.smooth(maxima, nPasses);
		updateDisplayRange();
	}

	public float getMaximaThreshold() {
		return maximaThreshold;
	}

	public void setMaximaThreshold(float maximaThreshold) {
		this.maximaThreshold = maximaThreshold;
	}

	public float getDisplayedMinimum() {
		return displayedMinimum;
	}

	public float getDisplayedMaximum() {
		return displayedMaximum;
	}

	public float getCurrentMinimum() {
		float min = maxima[0] & 0xffff;
		for(int i = 1; i < maxima.length; i++) {
			float v = maxima[i] & 0xffff;
			if(v < min)
				min = v;
		}
		return min;
	}

	public float getCurrentMaximum() {
		float max = maxima[0] & 0xffff;
		for(int i = 1; i < maxima.length; i++) {
			float v = maxima[i] & 0xffff;
			if(v > max)
				max = maxima[i];
		}
		return max;
	}

	public void setDisplayedMinimum(float min) {
		displayedMinimum = min;
		updateDisplayRange();
	}

	public void setDisplayedMaximum(float max) {
		displayedMaximum = max;
		updateDisplayRange();
	}

	private void updateDisplayRange() {
		boolean[] isMax = smp.isMaximum(maxima);
		int[] overlay = null;
		if(showColorOverlay) {
			File cf = files[currentIdx];
			File colorfile = new File(new File(cf.getParentFile(), "contributions"), cf.getName());
			if(!colorfile.exists())
				colorfile = new File(cf.getParentFile(), "contributions.vertices");
			if(colorfile.exists()) {
				try {
					System.out.println("Loading " + colorfile.getAbsolutePath());
					overlay = SphericalMaxProjection.loadIntData(colorfile.getAbsolutePath(), maxima.length);
				} catch(IOException e) {
					e.printStackTrace();
					showColorOverlay = false;
				}
			}
		}

		Point3f[] newvertices = new Point3f[smp.getSphere().nVertices];
		for(int i = 0; i < mesh.colors.length; i++) {
			float m = maxima[i] & 0xffff;
			if(showMaxima && isMax[i] && m > maximaThreshold)
				mesh.colors[i].set(1, 0, 0);
			else {
				m = (m - displayedMinimum) / (displayedMaximum - displayedMinimum);
				if(m < 0) m = 0;
				if(m > 1) m = 1;
				newvertices[i] = new Point3f(
					smp.getCenter().x + (1 + elevationFactor * m) * (mesh.vertices[i].x - smp.getCenter().x),
					smp.getCenter().y + (1 + elevationFactor * m) * (mesh.vertices[i].y - smp.getCenter().y),
					smp.getCenter().z + (1 + elevationFactor * m) * (mesh.vertices[i].z - smp.getCenter().z));
				if(showColorOverlay) {
					int rgb = overlay[i];
					int r = (rgb & 0xff0000) >> 16;
					int g = (rgb & 0xff00) >> 8;
					int b = (rgb & 0xff);

					double v = 255 * m;
					double c = 1 - (1/2.0 + m/2.0);
					r = Math.min(255 , (int)Math.round((c * r) + (1-c) * v));
					g = Math.min(255 , (int)Math.round((c * g) + (1-c) * v));
					b = Math.min(255 , (int)Math.round((c * b) + (1-c) * v));

					mesh.colors[i].set(r / 255f, g / 255f, b / 255f);
				} else if(lut == null) {
					mesh.colors[i].set(m, m, m);
				} else {
					int idx = Math.round(m * 255);
					int rv = lut[0][idx] & 0xff;
					int gv = lut[1][idx] & 0xff;
					int bv = lut[2][idx] & 0xff;
					mesh.colors[i].set(rv / 255f, gv / 255f, bv / 255f);
				}
			}
		}
		((IndexedTriangleArray)mesh.getGeometry()).setCoordinates(0, newvertices);
		setColors(mesh.colors);
	}

	void readColors(File file) throws IOException {
		this.maxima = smp.loadMaxima(file.getAbsolutePath());
	}

	public void setColors(Color3f[] colors) {
		((IndexedTriangleArray)mesh.getGeometry()).setColors(0, colors);
	}

	public void scaleForAngle(float angleFactor) throws IOException {
		Point3f[] vertices = smp.getSphere().getVertices();
		float cz = smp.getCenter().z;
		float radius = smp.getRadius();
		float bg = SphericalMaxProjection.getMode(maxima);
		System.out.println("bg = " + bg);
		for(int i = 0; i < vertices.length; i++) {
			float m = maxima[i] & 0xffff;
			float alpha = (float)Math.acos(Math.abs(vertices[i].z - cz) / radius);
			float factor = (float)(1.0 / (Math.cos(angleFactor * alpha)));
			float toSubtract = maxima[i] - bg - 10 > 0 ? bg + 10 : maxima[i];
			maxima[i] = (short)((m - toSubtract) * factor + toSubtract);
		}
		updateDisplayRange();
	}

	// timeline stuff
	@Override public void addInstant(ContentInstant ci) {}
	@Override public void removeInstant(int t) {}
	@Override public void setShowAllTimepoints(boolean b) {}

	@Override public int getNumberOfInstants() {
		return files.length;
	}

	@Override public ContentInstant getInstant(int t) {
		try {
			readColors(files[t]);
			updateDisplayRange();
		} catch(Exception e) {
			throw new RuntimeException("Cannot load " + files[t], e);
		}
		return super.getInstant(0);
	}

	@Override public ContentInstant getCurrent() {
		return super.getInstant(0);
	}

	@Override public void showTimepoint(int t) {
		try {
			currentIdx = t;
			readColors(files[t]);
			updateDisplayRange();
		} catch(Exception e) {
			throw new RuntimeException("Cannot load " + files[t], e);
		}
	}

	@Override public void showTimepoint(int t, boolean force) {
		if(force)
			super.showTimepoint(t, true);
		showTimepoint(t);
	}

	@Override public boolean isVisibleAt(int t) {
		return t >= 0 && t < getNumberOfInstants();
	}

	@Override public int getStartTime() {
		return 0;
	}

	@Override public int getEndTime() {
		return files.length - 1;
	}


	protected static final class CustomIndexedTriangleMesh extends CustomTriangleMesh {

		private final Point3f[] vertices;
		private final Color3f[] colors;
		private int[] faces;

		CustomIndexedTriangleMesh(Point3f[] vertices, Color3f[] colors, int[] faces) {
			super(null);
			this.vertices = vertices;
			this.colors = colors;
			this.faces = faces;
			update();
		}

		CustomIndexedTriangleMesh createNormalizedVersion(Point3f center, float radius) {
			Point3f[] nvertices = new Point3f[vertices.length];
			for(int i = 0; i < nvertices.length; i++)
				nvertices[i] = new Point3f();
			((IndexedTriangleArray)getGeometry()).getCoordinates(0, nvertices);

			for(int i = 0; i < vertices.length; i++) {
				nvertices[i].set(
					(nvertices[i].x - center.x) / radius,
					(nvertices[i].y - center.y) / radius,
					(nvertices[i].z - center.z) / radius);
			}
			return new CustomIndexedTriangleMesh(nvertices, colors, faces);
		}

		public void exportToPLY(String path) throws IOException {
			PrintStream out = new PrintStream(new FileOutputStream(path));
			out.println("ply");
			out.println("format ascii 1.0");
			out.println("element vertex " + vertices.length);
			out.println("property float x");
			out.println("property float y");
			out.println("property float z");
			out.println("property uchar red");
			out.println("property uchar green");
			out.println("property uchar blue");
			out.println("element face " + (faces.length / 3));
			out.println("property list uchar int vertex_index");
			out.println("end_header");
			for(int i = 0; i < vertices.length; i++) {
				Point3f v = vertices[i];
				Color3f c = colors[i];
				int r = Math.round(255 * c.x);
				int g = Math.round(255 * c.y);
				int b = Math.round(255 * c.z);
				out.println(v.x + " " + v.y + " " + v.z + " " + r + " " + g + " " + b);
			}
			int fl = faces.length / 3;
			for(int i = 0; i < fl; i++) {
				int f1 = faces[3 * i + 0];
				int f2 = faces[3 * i + 1];
				int f3 = faces[3 * i + 2];
				out.println("3 " + f1 + " " + f2 + " " + f3);
			}
			out.close();
		}

		@Override
		protected GeometryArray createGeometry() {
			if(vertices != null)
				return createGeometry(vertices, colors, faces);
			return null;
		}

		@Override
		public void calculateMinMaxCenterPoint(Point3f min,
					Point3f max, Point3f center) {

			min.x = min.y = min.z = Float.MAX_VALUE;
			max.x = max.y = max.z = Float.MIN_VALUE;
			for(int i = 0; i < faces.length; i++) {
				Point3f p = vertices[faces[i]];
				if(p.x < min.x) min.x = p.x;
				if(p.y < min.y) min.y = p.y;
				if(p.z < min.z) min.z = p.z;
				if(p.x > max.x) max.x = p.x;
				if(p.y > max.y) max.y = p.y;
				if(p.z > max.z) max.z = p.z;
			}
			center.x = (max.x + min.x) / 2;
			center.y = (max.y + min.y) / 2;
			center.z = (max.z + min.z) / 2;
		}

		protected GeometryArray createGeometry(Point3f[] vertices, Color3f[] colors, int[] faces) {

			IndexedTriangleArray ta = new IndexedTriangleArray(
				vertices.length,
					TriangleArray.COORDINATES |
					TriangleArray.COLOR_3, // |
					// TriangleArray.NORMALS,
				faces.length);

			ta.setCoordinates(0, vertices);
			ta.setColors(0, colors);

			ta.setCoordinateIndices(0, faces);
			ta.setColorIndices(0, faces);

			ta.setCapability(GeometryArray.ALLOW_COORDINATE_WRITE);
			ta.setCapability(GeometryArray.ALLOW_COLOR_WRITE);
			ta.setCapability(GeometryArray.ALLOW_INTERSECT);

			return ta;
		}

		private final Point2d p2d = new Point2d();
		private boolean roiContains(Point3f p, Transform3D volToIP, Canvas3D canvas, Polygon polygon) {
			Point3d locInImagePlate = new Point3d(p);
			volToIP.transform(locInImagePlate);
			canvas.getPixelLocationFromImagePlate(locInImagePlate, p2d);
			return polygon.contains(p2d.x, p2d.y);
		}

		@Override
		public void retain(Canvas3D canvas, Polygon polygon) {
			Transform3D volToIP = new Transform3D();
			canvas.getImagePlateToVworld(volToIP);
			volToIP.invert();

			Transform3D toVWorld = new Transform3D();
			this.getLocalToVworld(toVWorld);
			volToIP.mul(toVWorld);

			ArrayList<Integer> f = new ArrayList<Integer>();
			for(int i = 0; i < faces.length; i += 3) {
				Point3f p1 = vertices[faces[i]];
				Point3f p2 = vertices[faces[i + 1]];
				Point3f p3 = vertices[faces[i + 2]];
				if(roiContains(p1, volToIP, canvas, polygon) ||
						roiContains(p2, volToIP, canvas, polygon) ||
						roiContains(p3, volToIP, canvas, polygon)) {
					f.add(faces[i]);
					f.add(faces[i + 1]);
					f.add(faces[i + 2]);
				}
			}
			faces = new int[f.size()];
			for(int i = 0; i < faces.length; i++)
				faces[i] = f.get(i);
			update();
		}
	}
}

