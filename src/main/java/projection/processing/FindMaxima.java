package projection.processing;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import javax.vecmath.Point2f;
import javax.vecmath.Point3f;

public class FindMaxima {

	/**
	 * the following constants are used to set bits corresponding to pixel
	 * types
	 */
	private final static byte MAXIMUM    = (byte)  1; // marks local maxima (irrespective of noise tolerance)
	private final static byte LISTED     = (byte)  2; // marks points currently in the list
	private final static byte PROCESSED  = (byte)  4; // marks points processed previously
	private final static byte MAX_AREA   = (byte)  8; // marks areas near a maximum, within the tolerance
	private final static byte EQUAL      = (byte) 16; // marks contigous maximum points of equal level
	private final static byte MAX_POINT  = (byte) 32; // marks a single point standing for a maximum

	private final SphericalMaxProjection smp;
	private final int[][] neighbors;

	public FindMaxima(SphericalMaxProjection smp) {
		this.smp = smp;
		int[] faces = smp.getSphere().faces;
		this.neighbors = calculateNeighbors(faces, smp.getSphere().nVertices);
	}

	@SuppressWarnings("unchecked")
	private int[][] calculateNeighbors(int[] faces, int nVertices) {
		HashSet<Integer>[] set = new HashSet[nVertices];
		for(int i = 0; i < nVertices; i++)
			set[i] = new HashSet<Integer>();

		for(int i = 0; i < faces.length; i += 3) {
			int f1 = faces[i];
			int f2 = faces[i + 1];
			int f3 = faces[i + 2];
			set[f1].add(f2);
			set[f1].add(f3);

			set[f2].add(f1);
			set[f2].add(f3);

			set[f3].add(f1);
			set[f3].add(f2);
		}

		int[][] neigh = new int[nVertices][];
		for(int i = 0; i < nVertices; i++) {
			neigh[i] = new int[set[i].size()];
			int n = 0;
			for(int neighbor : set[i])
				neigh[i][n++] = neighbor;
		}
		return neigh;
	}

	public ArrayList<Point3f> findMaxima(short[] values, boolean[] isVertexMax, double tolerance) {
		byte[] types = new byte[values.length];
		long[] maxPoints = getSortedMaxPoints(values, types);
		return analyzeAndMarkMaxima(values, types, maxPoints, isVertexMax, tolerance);
	}


	long[] getSortedMaxPoints(short[] values, byte[] types) {
		int nMax = 0; // counts local maxima
		for(int i = 0; i < values.length; i++) {
			int v = (values[i] & 0xffff);
			boolean isMax = true;
			for (int n : neighbors[i]) {
				int vNeighbor = values[n] & 0xffff;
				if (vNeighbor > v) {
					isMax = false;
					break;
				}
			}
			if (isMax) {
				types[i] = MAXIMUM;
				nMax++;
			}
		}

		// upper 32 bit is pixel value,
		// lower 32 bit is pixel pos
		long[] maxPoints = new long[nMax];
		int iMax = 0;
		for(int i = 0; i < values.length; i++) {
			if (types[i] == MAXIMUM) {
				int iValue = values[i] & 0xffff;
				maxPoints[iMax++] = (long) iValue << 32	| i;
			}
		}
		Arrays.sort(maxPoints);
		return maxPoints;
	}

	ArrayList<Point3f> analyzeAndMarkMaxima(
			short[] values,
			byte[] types,
			long[] maxPoints,
			boolean[] isVertexMax,
			double tolerance) {

		int nMax = maxPoints.length;
		int[] pList = new int[values.length];   // here we enter points
							// starting from a
							// maximum
		ArrayList<Point3f> maxima = new ArrayList<Point3f>();

		Point3f[] vertices = smp.getSphere().vertices;
		Point2f tmp = new Point2f();
		// process all maxima, starting with the highest
		for (int iMax = nMax - 1; iMax >= 0; iMax--) {
			int offset0 = (int) maxPoints[iMax];
			if ((types[offset0] & PROCESSED) != 0)
				continue;
			// we create a list of connected points and start the
			// list at the current maximum
			int v0 = values[offset0] & 0xffff;
			pList[0] = offset0;
			// mark first point as equal height (to itself) and listed
			types[offset0] |= (EQUAL | LISTED);
			int listLen = 1; // number of elements in the list
			int listI = 0; // index of current element in the list
			boolean maxPossible = true; // it may be a true maximum
			smp.getPolar(vertices[offset0], tmp);
			double lonEqual = tmp.x;
			double latEqual = tmp.y;
			int nEqual = 1; // counts xEqual/yEqual points that we
					// use for averaging

			do { // while neigbor list is not fully processed (to listLen)

				int offset = pList[listI];
				// analyze all neighbors at the same level
				for (int offset2 : neighbors[offset]) {
					if ((types[offset2] & LISTED) == 0) {
						if ((types[offset2] & PROCESSED) != 0) {
							// reached a previously processed point:
							// -> no maximum
							maxPossible = false;
							break;
						}
						int v2 = values[offset2] & 0xffff;
						if (v2 > v0) {
							// reached a higher point: -> no maximum
							maxPossible = false;
							break;
						} else if (v2 >= v0 - (float) tolerance) {
							pList[listLen] = offset2;
							listLen++;
							// we have found a new point
							// within the tolerance
							types[offset2] |= LISTED;
							if (v2 == v0) {
								// prepare finding center of equal points
								// (in case single point needed)
								types[offset2] |= EQUAL;
								smp.getPolar(vertices[offset2], tmp);
								lonEqual += tmp.x;
								latEqual += tmp.y;
								nEqual++;
							}
						}
					} // if isWithin & not LISTED
				} // for directions d
				listI++;
			} while (listI < listLen);

			int resetMask = ~(maxPossible ? LISTED
					: (LISTED | EQUAL));
			lonEqual /= nEqual;
			latEqual /= nEqual;
			double minDist2 = 1e20;
			int nearestI = 0;
			for (listI = 0; listI < listLen; listI++) {
				int offset = pList[listI];
				smp.getPolar(vertices[offset], tmp);
				double lon = tmp.x;
				double lat = tmp.y;
				types[offset] &= resetMask; // reset attributes no longer needed
				types[offset] |= PROCESSED; // mark as processed
				if (maxPossible) {
					types[offset] |= MAX_AREA;
					if ((types[offset] & EQUAL) != 0) {
						double dist2 = (lonEqual - lon)
								* (lonEqual - lon)
								+ (latEqual - lat)
								* (latEqual - lat);
						if (dist2 < minDist2) {
							minDist2 = dist2;
							// this could be the best
							// "single maximum" point
							nearestI = listI;
						}
					}
				}
			} // for listI
			if (maxPossible) {
				int offset = pList[nearestI];
				types[offset] |= MAX_POINT;
				maxima.add(new Point3f(vertices[offset]));
				if(isVertexMax != null)
					isVertexMax[offset] = true;
			}
		} // for all maxima iMax
		return maxima;
	} // void analyzeAndMarkMaxima
}
