package projection.util;

import ij.ImagePlus;
import ij.measure.Calibration;
import ij.process.ImageProcessor;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Stack;

import javax.vecmath.Point3f;
import javax.vecmath.Point3i;
import javax.vecmath.Vector3f;

public class MeshToMask {

	private static final int FG = 255;

	private final Image image;

	public MeshToMask(ImagePlus image) {
		this.image = new Image(image, false, true);
	}

	public void drawTriangle(Point3f a, Point3f b, Point3f c, int v) {
		image.drawTriangle(a, b, c, v);
	}

	public void drawLine(Point3f a, Point3f b, int v) {
		image.drawLine(a, b, v);
	}

	public void drawMesh(List<Point3f> mesh, int v) {
		image.drawMesh(mesh, v);
	}

	public static List<Point3i> collectSurfacePoints(
			ImagePlus into, List<Point3f> mesh) {
		Image im = new Image(into, true, false);
		im.drawMesh(mesh, FG);
		return new ArrayList<Point3i>(im.points);
	}

	public static List<Point3i> collectVertexPoints(
			ImagePlus into, List<Point3f> mesh) {
		Set<Point3f> vertices = new HashSet<Point3f>(mesh);
		List<Point3i> list = new ArrayList<Point3i>(vertices.size());
		Calibration c = into.getCalibration();
		for(Point3f p : vertices)
			list.add(new Point3i(
				(int)Math.round(p.x / c.pixelWidth),
				(int)Math.round(p.y / c.pixelHeight),
				(int)Math.round(p.z / c.pixelDepth)));
		return list;
	}

	public static List<Point3i> collectEdgePoints(
			ImagePlus into, List<Point3f> mesh) {

		Image im = new Image(into, true, false);
		for(int i = 0; i < mesh.size(); i += 3) {
			Point3f p1 = mesh.get(i);
			Point3f p2 = mesh.get(i);
			Point3f p3 = mesh.get(i);
			im.drawLine(p1, p2, FG);
			im.drawLine(p2, p3, FG);
			im.drawLine(p3, p1, FG);
		}
		return new ArrayList<Point3i>(im.points);
	}

	public static void drawMesh(ImagePlus into, List<Point3f> mesh) {
		drawMesh(into, mesh, FG);
	}

	public static void drawMesh(ImagePlus into, List<Point3f> mesh, int v) {
		Image im = new Image(into, false, true);
		im.drawMesh(mesh, v);
	}

	public static void fillMesh(ImagePlus into, List<Point3f> mesh) {
		fillMesh(into, mesh, FG);
	}

	public static void drawMeshVertices(ImagePlus into, List<Point3f> mesh) {
		drawMeshVertices(into, mesh, FG);
	}

	public static void drawMeshVertices(ImagePlus into, List<Point3f> mesh, int v) {
		Image im = new Image(into, false, true);
		im.drawMeshVertices(mesh, v);
	}

	public static void fillMesh(ImagePlus into, List<Point3f> mesh, int v) {
		// try to guess a point inside the mesh
		Point3f p1 = mesh.get(0);
		Point3f p2 = mesh.get(1);
		Point3f p3 = mesh.get(2);
		Vector3f v1 = new Vector3f();
		Vector3f v2 = new Vector3f();
		v1.sub(p2, p1);
		v2.sub(p3, p1);
		v1.cross(v1, v2);
		v1.normalize();

		Image im = new Image(into, false, true);
		v1.x *= (3 * im.pw);
		v1.y *= (3 * im.ph);
		v1.z *= (3 * im.pd);

		Point3f p = new Point3f(p1);
		p.add(p2);
		p.add(p3);
		p.scale(1.0f / 3);
		p.sub(v1);

		fillMesh(into, mesh, p, v);
	}

	public static void fillMesh(ImagePlus into,
			List<Point3f> mesh, Point3f src, int v) {

		Image im = new Image(into, false, true);
		im.drawMesh(mesh, v);
		im.fillMesh(src, v);
	}

	private static final class Image {

		private final int w, h, d;
		private final double pw, ph, pd;
		private final double mininc;
		private final ImageProcessor[] data;
		private final boolean collect;
		private final boolean draw;
		private final Set<Point3i> points;

		public Image(ImagePlus imp, boolean collect, boolean draw) {
			this.draw = draw;
			this.collect = collect;

			points = collect ? new HashSet<Point3i>() : null;

			this.w = imp.getWidth();
			this.h = imp.getHeight();
			this.d = imp.getStackSize();
			this.data = new ImageProcessor[d];
			for(int z = 0; z < d; z++)
				data[z] = imp.getStack().getProcessor(z + 1);

			Calibration cal = imp.getCalibration();
			pw = cal.pixelWidth;
			ph = cal.pixelHeight;
			pd = cal.pixelDepth;

			mininc = Math.min(pw, Math.min(ph, pd)) * 0.5;
		}

		void set(int x, int y, int z, int v) {
			if(x < 0 || x >= w || y < 0 || y >= h || z < 0 || z >= d)
				return;
			if(draw) {
				int i = y * w + x;
				data[z].set(i, v);
			}

			if(collect)
				points.add(new Point3i(x, y, z));
		}

		int get(int x, int y, int z) {
			return data[z].get(x, y);
		}

		void drawLine(Point3f p1, Point3f p2, int v) {

			Vector3f n = new Vector3f();
			n.sub(p1, p2);
			double d = n.length();
			n.normalize();

			Point3f p = new Point3f();
			for (double i = 0; i <= d; i += mininc) {
				p.scaleAdd((float)i, n, p2);
				set((int) Math.floor(p.x / pw +.5),
						(int) Math.floor(p.y / ph +.5),
						(int) Math.floor(p.z / pd +.5),
						v);
			}
		}

		void drawMeshVertices(List<Point3f> mesh, int v) {
			for(Point3f p : mesh)
				set((int)Math.round(p.x / pw),
					(int)Math.round(p.y / ph),
					(int)Math.round(p.z / pd),
					v);
		}

		private final Vector3f nTmp = new Vector3f();
		private final Point3f pTmp = new Point3f();

		void drawTriangle(Point3f p1, Point3f p2, Point3f p3, int v) {
			nTmp.sub(p1, p2);
			double d = nTmp.length();
			nTmp.normalize();

			for (double j = 0; j <= d; j += mininc) {
				pTmp.scaleAdd((float)j, nTmp, p2);
				drawLine(pTmp, p3, v);
			}
		}

		void drawMesh(List<Point3f> mesh, int v) {
			for(int i = 0; i < mesh.size(); i += 3) {
				Point3f p1 = mesh.get(i);
				Point3f p2 = mesh.get(i + 1);
				Point3f p3 = mesh.get(i + 2);

				drawTriangle(p1, p2, p3, v);
			}
		}

		void fillMesh(Point3f source, int v) {
			fillMesh(new Point3i(
				(int)Math.round(source.x / pw),
				(int)Math.round(source.y / ph),
				(int)Math.round(source.z / pd)), v);
		}

		void fillMesh(Point3i source, int v) {
			Stack<Point3i> s= new Stack<Point3i>();
			s.push(source);
			while(!s.isEmpty()) {
				Point3i p = s.pop();

				if(get(p.x, p.y, p.z) == v)
					continue;

				set(p.x, p.y, p.z, v);

				if(p.x > 0)
					s.push(new Point3i(p.x - 1, p.y, p.z));
				if(p.x < w - 1)
					s.push(new Point3i(p.x + 1, p.y, p.z));
				if(p.y > 0)
					s.push(new Point3i(p.x, p.y - 1, p.z));
				if(p.y < h - 1)
					s.push(new Point3i(p.x, p.y + 1, p.z));
				if(p.z > 0)
					s.push(new Point3i(p.x, p.y, p.z - 1));
				if(p.z < d - 1)
					s.push(new Point3i(p.x, p.y, p.z + 1));
			}
		}
	}
}

