package projection.processing;

import ij.ImagePlus;
import ij.ImageStack;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.scijava.vecmath.Point3f;
import org.scijava.vecmath.Tuple3f;
import org.scijava.vecmath.Vector3f;

import projection.util.MeshToMask;


public class IndexedTriangleMesh {

	protected int[] faces;
	protected Point3f[] vertices;

	public final int nVertices;
	public final int nFaces;

	private IndexedTriangleMesh(int nFaces, int nVertices) {
		this.nFaces = nFaces;
		this.nVertices = nVertices;
	}

	public IndexedTriangleMesh(Point3f[] vertices, int[] faces) {
		this.nFaces = faces.length;
		this.nVertices = vertices.length;
		this.faces = faces;
		this.vertices = vertices;
	}

	public IndexedTriangleMesh(List<Point3f> mesh) {
		Map<Point3f, Integer> vertexToIndex =
				new HashMap<Point3f, Integer>();

		this.nFaces = mesh.size();
		this.faces = new int[nFaces];
		List<Point3f> v = new ArrayList<Point3f>();

		for(int i = 0; i < nFaces; i++) {
			Point3f p = mesh.get(i);

			if(!vertexToIndex.containsKey(p)) {
				Point3f newp = new Point3f(p);
				vertexToIndex.put(newp, v.size());
				v.add(newp);
			}
			faces[i] = vertexToIndex.get(p);
		}
		this.nVertices = v.size();
		this.vertices = new Point3f[nVertices];
		v.toArray(vertices);
	}

	public int[] getFaces() {
		return faces;
	}

	public Point3f[] getVertices() {
		return vertices;
	}

	public List<Point3f> createMesh() {
		List<Point3f> mesh = new ArrayList<Point3f>(faces.length);
		for(int i = 0; i < faces.length; i++)
			mesh.add(vertices[faces[i]]);
		return mesh;
	}

	public Vector3f[] createNormals() {
		Vector3f[] normals = new Vector3f[nVertices];
		for(int i = 0; i < nVertices; i++)
			normals[i] = new Vector3f();

		Vector3f v1 = new Vector3f(), v2 = new Vector3f();
		for(int i = 0; i < nFaces; i += 3) {
			int f1 = faces[i];
			int f2 = faces[i + 1];
			int f3 = faces[i + 2];

			v1.sub(vertices[f2], vertices[f1]);
			v2.sub(vertices[f3], vertices[f1]);
			v1.cross(v1, v2);

			normals[f1].add(v1);
			normals[f2].add(v1);
			normals[f3].add(v1);
		}
		for(Vector3f v : normals)
			v.normalize();
		return normals;
	}

	public void getCenterOfGravity(Tuple3f cog) {
		cog.set(0, 0, 0);
		for(Point3f v : vertices)
			cog.add(v);
		cog.scale(1f / vertices.length);
	}

	public ImagePlus createOverlay(ImagePlus image) {
		return createOverlay(image, 0x00ff00);
	}

	public ImagePlus createOverlay(ImagePlus image, int col) {
		ImageStack stack = new ImageStack(
			image.getWidth(), image.getHeight());
		int d = image.getStackSize();
		for(int z = 0; z < d; z++)
			stack.addSlice("", image.getStack().
				getProcessor(z + 1).convertToRGB());
		ImagePlus result = new ImagePlus("result", stack);
		result.setCalibration(image.getCalibration());

		MeshToMask.drawMesh(result, createMesh(), col);
		return result;
	}

	@Override
	public Object clone() {
		IndexedTriangleMesh cp = new IndexedTriangleMesh(
				nFaces, nVertices);
		cp.faces = new int[nFaces];
		System.arraycopy(faces, 0, cp.faces, 0, nFaces);

		cp.vertices = new Point3f[nVertices];
		for(int i = 0; i < nVertices; i++)
			cp.vertices[i] = new Point3f(vertices[i]);
		return cp;
	}
}

