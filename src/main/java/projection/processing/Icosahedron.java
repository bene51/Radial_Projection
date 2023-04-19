package projection.processing;

import ij3d.Image3DUniverse;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.scijava.vecmath.Color3f;
import org.scijava.vecmath.Point2i;
import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Vector3f;

import customnode.CustomPointMesh;
import customnode.CustomTriangleMesh;

public class Icosahedron extends IndexedTriangleMesh {

	private static final float TAO = (float)((1 + Math.sqrt(5)) / 2.0);

	public static Point3f[] createVertices(float r) {
		Point3f[] v = new Point3f[] {
			new Point3f(1, TAO, 0), new Point3f(-1, TAO, 0), new Point3f(1, -TAO, 0), new Point3f(-1, -TAO, 0),
			new Point3f(0, 1, TAO), new Point3f(0, -1, TAO), new Point3f(0, 1, -TAO), new Point3f(0, -1, -TAO),
			new Point3f(TAO, 0, 1), new Point3f(-TAO, 0, 1), new Point3f(TAO, 0, -1), new Point3f(-TAO, 0, -1)};
		float scale = (float)Math.sqrt((r * r) / (1f + TAO * TAO));
		for(Point3f p : v)
			p.scale(scale);
		return v;
	}

	public static int[] createFaces() {
		return new int[] {
			0, 1, 4, 1, 9, 4, 4, 9, 5, 5, 9, 3,
			2, 3, 7, 3, 2, 5, 7, 10, 2, 0, 8, 10,
			0, 4, 8, 8, 2, 10, 8, 4, 5, 8, 5, 2,
			1, 0, 6, 11, 1, 6, 3, 9, 11, 6, 10, 7,
			3, 11, 7, 11, 6, 7, 6, 0, 10, 9, 1, 11};
	}

	public Icosahedron(float r) {
		super(createVertices(r), createFaces());
	}

	public IndexedTriangleMesh createFlatVersion(float s) {
		float h = (float)(0.5 * Math.sqrt(3) * s);
		Point3f[] vertices = new Point3f[22];
		vertices[0] = new Point3f((1 + 0) * s, 0, 0);
		vertices[1] = new Point3f((1 + 1) * s, 0, 0);
		vertices[2] = new Point3f((1 + 2) * s, 0, 0);
		vertices[3] = new Point3f((1 + 3) * s, 0, 0);
		vertices[4] = new Point3f((1 + 4) * s, 0, 0);

		vertices[5]  = new Point3f((.5f + 0) * s, h, 0);
		vertices[6]  = new Point3f((.5f + 1) * s, h, 0);
		vertices[7]  = new Point3f((.5f + 2) * s, h, 0);
		vertices[8]  = new Point3f((.5f + 3) * s, h, 0);
		vertices[9]  = new Point3f((.5f + 4) * s, h, 0);
		vertices[10] = new Point3f((.5f + 5) * s, h, 0);

		vertices[11] = new Point3f(0 * s, 2 * h, 0);
		vertices[12] = new Point3f(1 * s, 2 * h, 0);
		vertices[13] = new Point3f(2 * s, 2 * h, 0);
		vertices[14] = new Point3f(3 * s, 2 * h, 0);
		vertices[15] = new Point3f(4 * s, 2 * h, 0);
		vertices[16] = new Point3f(5 * s, 2 * h, 0);

		vertices[17] = new Point3f((0.5f + 0) * s, 3 * h, 0);
		vertices[18] = new Point3f((0.5f + 1) * s, 3 * h, 0);
		vertices[19] = new Point3f((0.5f + 2) * s, 3 * h, 0);
		vertices[20] = new Point3f((0.5f + 3) * s, 3 * h, 0);
		vertices[21] = new Point3f((0.5f + 4) * s, 3 * h, 0);

		int[] faces = new int[] {
			5, 12, 6, 12, 13, 6, 6, 13, 7, 7, 13, 14,
			8, 14, 15, 14, 8, 7, 15, 9, 8, 10, 4, 9,
			5, 6, 0, 3, 8, 9, 1, 6, 7, 2, 7, 8,
			12, 5, 11, 17, 12, 11, 14, 13, 19, 16, 9, 15,
			14, 20, 15, 21, 16, 15, 16, 10, 9, 13, 12, 18};

		return new IndexedTriangleMesh(vertices, faces);
	}

	private static int getMiddlePoint(ArrayList<Point3f> vertices, Point2i v1v2, HashMap<Point2i, Integer> cache) {
		if(v1v2.x > v1v2.y)
			v1v2.set(v1v2.y, v1v2.x);

		if(cache.containsKey(v1v2))
			return cache.get(v1v2);

		Point3f p1 = vertices.get(v1v2.x);
		Point3f p2 = vertices.get(v1v2.y);

		Point3f middle = new Point3f(p1);
		middle.add(p2);
		middle.scale(0.5f);

		int i = vertices.size();
		vertices.add(middle);
		cache.put(new Point2i(v1v2), i);
		return i;
	}

	/**
	 * Based on http://blog.andreaskahler.com/2009/06/creating-icosphere-mesh-in-code.html
	 * @param r
	 * @param levels Number of levels to split each edge.
	 * @return
	 */
	public IndexedTriangleMesh createBuckyballFast(float r, int levels) {
		ArrayList<Point3f> vertices = new ArrayList<Point3f>(Arrays.asList(this.vertices));
		int[] faces = this.faces;

		Point2i v1v2 = new Point2i();
		for(int l = 0; l < levels; l++) {
			HashMap<Point2i, Integer> cache = new HashMap<Point2i, Integer>();
			// ArrayList<Integer> facesN = new ArrayList<Integer>();
			int[] facesN = new int[faces.length * 4];
			for(int i = 0; i < faces.length; i += 3) {
				int v1 = faces[i];
				int v2 = faces[i + 1];
				int v3 = faces[i + 2];

				v1v2.set(v1, v2); int a = getMiddlePoint(vertices, v1v2, cache);
				v1v2.set(v2, v3); int b = getMiddlePoint(vertices, v1v2, cache);
				v1v2.set(v3, v1); int c = getMiddlePoint(vertices, v1v2, cache);

				facesN[i * 4 +  0] = v1;
				facesN[i * 4 +  1] = a;
				facesN[i * 4 +  2] = c;
				facesN[i * 4 +  3] = v2;
				facesN[i * 4 +  4] = b;
				facesN[i * 4 +  5] = a;
				facesN[i * 4 +  6] = v3;
				facesN[i * 4 +  7] = c;
				facesN[i * 4 +  8] = b;
				facesN[i * 4 +  9] = a;
				facesN[i * 4 + 10] = b;
				facesN[i * 4 + 11] = c;
			}
			faces = facesN;
		}

		Point3f[] v = new Point3f[vertices.size()];
		vertices.toArray(v);

		Point3f m = new Point3f();
		Vector3f d = new Vector3f();
		for(Point3f vtx : v) {
			d.sub(vtx, m);
			d.normalize();
			vtx.scaleAdd(r, d, m);
		}
		return new IndexedTriangleMesh(v, faces);
	}

	public static void main(String[] args) {


		Icosahedron icosa = new Icosahedron(1f);
		// univ.addCustomMesh(new CustomTriangleMesh(icosa.createMesh()), "icosahedron");

		float r = 3f;
		int subd = 4;
		int levels = (int)Math.floor(Math.log(subd) / Math.log(2));
		System.out.println(levels + " levels");

		long start = System.currentTimeMillis();
		IndexedTriangleMesh bucky = icosa.createBuckyballFast(r, levels);
		long end = System.currentTimeMillis();
		System.out.println("Approach 2 took " + (end - start) + " ms (" + bucky.nVertices + " vtx)");

		Image3DUniverse univ = new Image3DUniverse();
		univ.show();

		List<Point3f> spherepoints = Arrays.asList(bucky.getVertices());
		List<Point3f> triangles = bucky.createMesh();

		CustomPointMesh cpm = new CustomPointMesh(spherepoints);
		cpm.setColor(new Color3f(1, 0, 0));
		cpm.setPointSize(3f);
		univ.addCustomMesh(cpm, "Sphere points");

		CustomTriangleMesh ctm = new CustomTriangleMesh(triangles);
		ctm.setColor(new Color3f(0.4f, 0.4f, 1));
		univ.addCustomMesh(ctm, "triangles");
	}
}
