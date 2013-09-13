package projection.util;

import javax.vecmath.Point3f;

public final class PointMatch implements Comparable<PointMatch> {
	public final Point3f p1;
	public final Point3f p2;
	public final float distance2;

	public PointMatch(Point3f p1, Point3f p2) {
		this.p1 = p1;
		this.p2 = p2;
		this.distance2 = p1.distanceSquared(p2);
	}

	public boolean equals(Object o) {
		PointMatch pm = (PointMatch)o;
		return pm.p1.equals(p1) && pm.p2.equals(p2);
	}

	public int hashCode() {
		return p1.hashCode() * p2.hashCode();
	}

	public int compareTo(PointMatch o) {
		if(this.distance2 < o.distance2)
			return -1;
		if(this.distance2 > o.distance2)
			return +1;
		return 0;
	}
}

