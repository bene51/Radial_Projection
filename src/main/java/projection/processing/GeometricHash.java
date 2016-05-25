package projection.processing;

import java.util.ArrayList;
import java.util.Arrays;

import javax.vecmath.Matrix4f;
import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

import projection.util.PointMatch;
import fiji.util.KDTree;
import fiji.util.NNearestNeighborSearch;
import fiji.util.node.Leaf;

public class GeometricHash {

	public static void main(String[] args) {
		ArrayList<Point3f> pts = new ArrayList<Point3f>();
		for(int i = 0; i < 4; i++) {
			pts.add(new Point3f(10 * (float)Math.random(), 10 * (float)Math.random(), 10 * (float)Math.random()));
		}
		System.out.println("pts: " + pts);
		hash(pts, false);
	}

	public static ArrayList<PointMatch> establishCorrespondences(ArrayList<Point3f> pts1, ArrayList<Point3f> pts2) {
		ArrayList<PointMatch> candidates = new ArrayList<PointMatch>();
		if(pts1.size() < 4 || pts2.size() < 4)
			return candidates;
		ArrayList<NodeND> params1 = hash(pts1, false);
		ArrayList<NodeND> params2 = hash(pts2, false);

		KDTree<NodeND> tree = new KDTree<NodeND>(params1);
		NNearestNeighborSearch<NodeND> nn = new NNearestNeighborSearch<NodeND>(tree);

		for(NodeND p : params2) {
			NodeND[] neighbors = nn.findNNearestNeighbors(p, 2);
			NodeND no1 = neighbors[0];
			NodeND no2 = neighbors[1];
			double distNo1 = p.distanceTo(no1);
			double distNo2 = p.distanceTo(no2);
			System.out.println("distNo1 = " + distNo1 + " distNo2 = " + distNo2);
			if(2 * distNo1 <= distNo2)
				candidates.add(new PointMatch(p.point, no1.point));
		}

		return candidates;
	}

	public static ArrayList<NodeND> hash(ArrayList<Point3f> pts, boolean includePosition) {
		ArrayList<Node3D> nodes = new ArrayList<Node3D>(pts.size());
		for(int i = 0; i < pts.size(); i++)
			nodes.add(new Node3D(pts.get(i), i));
		KDTree<Node3D> tree = new KDTree<Node3D>(nodes);
		NNearestNeighborSearch<Node3D> nn = new NNearestNeighborSearch<Node3D>(tree);

		ArrayList<NodeND> ret = new ArrayList<NodeND>(pts.size());
		int dim = includePosition ? 9 : 6;
		for(int i = 0; i < nodes.size(); i++) {
			int pIdx = 0;
			Node3D a = nodes.get(i);
			float[] param = new float[dim];
			Node3D[] neighbors = nn.findNNearestNeighbors(a, 4);
			Node3D b = neighbors[1];
			Node3D c = neighbors[2];
			Node3D d = neighbors[3];

			if(includePosition) {
				param[pIdx++] = a.p.x;
				param[pIdx++] = a.p.y;
				param[pIdx++] = a.p.z;
			}

			// the normalized vector ad is the local x-axis
			Vector3f x = new Vector3f(d.p);
			x.sub(a.p);
			x.normalize();
			param[pIdx++] = a.p.distance(d.p);

			// n is the normalized vector ad, which is the normal vector
			// of a plane perpendicular to ad, i.e. the local x vector. In this
			// plane is the local y- and z-axis.
			Vector3f ab = new Vector3f(b.p);
			ab.sub(a.p);
			Vector3f n = new Vector3f(x);
			float distB = n.dot(ab);
			n.scale(distB);
			Point3f bDash = new Point3f(b.p);
			bDash.sub(n);
			Vector3f y = new Vector3f(bDash);
			y.sub(a.p);
			y.normalize();
			param[pIdx++] = distB;
			param[pIdx++] = bDash.distance(a.p);

			// the local z axis is perpendicular to the other axis:
			Vector3f z = new Vector3f();
			z.cross(x, y);
			Point3f cDash = new Point3f(c.p);
			cDash.sub(a.p);
			Matrix4f m = new Matrix4f(
					x.x, y.x, z.x, 0,
					x.y, y.y, z.y, 0,
					x.z, y.z, z.z, 0,
					  0,   0,   0, 1);
			m.invert();
			m.transform(cDash);
			param[pIdx++] = cDash.x;
			param[pIdx++] = cDash.y;
			param[pIdx++] = cDash.z;
			ret.add(new NodeND(param, a.p, i));
		}

		return ret;
	}

	private static class Node3D implements Leaf<Node3D> {

		final Point3f p;
		final int index;

		public Node3D(final Point3f p, int index) {
			this.p = p;
			this.index = index;
		}

		@SuppressWarnings("unused")
		public Node3D(final Node3D node) {
			this.p = (Point3f)node.p.clone();
			this.index = node.index;
		}

		@Override
		public boolean isLeaf() {
			return true;
		}

		@SuppressWarnings("unused")
		public boolean equals(final Node3D o) {
	                 return p.equals(o.p);
		}

		@Override
		public float distanceTo(final Node3D o) {
			return p.distance(o.p);
		}

		@Override
		public float get(final int k) {
			switch(k) {
				case 0: return p.x;
				case 1: return p.y;
				case 2: return p.z;
			}
			return 0f;
		}

		@Override
		public String toString() {
			return index + ": " + p.toString();
		}

		@Override
		public Node3D[] createArray(final int n) {
			return new Node3D[n];
		}

		@Override
		public int getNumDimensions() {
			return 3;
		}
	}

	static class NodeND implements Leaf<NodeND> {

		final float[] param;
		final int index;
		final Point3f point;
		final int dimensions;

		public NodeND(float[] param, Point3f point, int index) {
			this.param = param;
			this.point = point;
			this.index = index;
			this.dimensions = param.length;
		}

		public NodeND(final NodeND node) {
			this.param = node.param.clone();
			this.point = node.point;
			this.index = node.index;
			this.dimensions = node.dimensions;
		}

		@Override
		public boolean isLeaf() {
			return true;
		}

		public boolean equals(final NodeND o) {
			return Arrays.equals(param, o.param);
		}

		@Override
		public float distanceTo(final NodeND o) {
			double dist = 0;
			for(int i = 0; i < dimensions; i++) {
				double d = o.param[i] - param[i];
				dist += d * d;
			}
			return (float)Math.sqrt(dist);
		}

		@Override
		public float get(final int k) {
			return param[k];
		}

		@Override
		public String toString() {
			return Arrays.toString(param);
		}

		@Override
		public NodeND[] createArray(final int n) {
			return new NodeND[n];
		}

		@Override
		public int getNumDimensions() {
			return dimensions;
		}
	}
}
