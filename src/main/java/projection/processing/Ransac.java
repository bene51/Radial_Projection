package projection.processing;

import java.util.ArrayList;
import java.util.HashSet;

import org.scijava.vecmath.Matrix4f;
import org.scijava.vecmath.Point3f;

import projection.util.PointMatch;

public class Ransac {

	private final ArrayList<PointMatch> candidates;
	private final Point3f cor;

	private final int iterations = 1000;
	private final float maxEpsilon = 50f;
	private final float minInlierRatio = 0.1f;

	public Ransac(ArrayList<PointMatch> candidates, Point3f cor) {
		this.candidates = candidates;
		this.cor = cor;
	}

	public Matrix4f run() {
		if(candidates.size() < 6)
			return null;

		Matrix4f bestRigid = new Matrix4f();
		int bestNumInlier = 0;
		boolean foundSomething = false;

		Matrix4f matrix = new Matrix4f();
		for(int it = 0; it < iterations; it++) {
			HashSet<Integer> inlierIndices = getNRandomIntegers(0, candidates.size(), 6);
			ArrayList<PointMatch> inliers = new ArrayList<PointMatch>();
			for(int inlierIndex : inlierIndices)
				inliers.add(candidates.get(inlierIndex));

			ICPRegistration.bestRigid(inliers, matrix, cor);

			Point3f tmp = new Point3f();
			int numInlier = 6;
			for(int i = 0; i < candidates.size(); i++) {
				if(inlierIndices.contains(i))
					continue;
				PointMatch pm = candidates.get(i);
				tmp.set(pm.p2);
				matrix.transform(tmp);
				double distance = tmp.distance(pm.p1);
				if(distance <= maxEpsilon)
					numInlier++;
			}
			if(numInlier >= minInlierRatio * candidates.size() && numInlier > bestNumInlier) {
				bestNumInlier = numInlier;
				bestRigid.set(matrix);
				foundSomething = true;
			}
			if(numInlier > bestNumInlier) {
				bestNumInlier = numInlier;
			}
		}
		System.out.println("just found maximally " + bestNumInlier + " inliers");
		return bestRigid;
	}

	private static final HashSet<Integer> getNRandomIntegers(int offs, int length, int n) {
		HashSet<Integer> found = new HashSet<Integer>();

		while(found.size() < n) {
			int c = offs + (int)Math.round((length - 1) * Math.random());
			found.add(c);
		}
		return found;
	}
}
