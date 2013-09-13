package projection.processing;

import ij.process.ColorProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.geom.GeneralPath;
import java.awt.geom.PathIterator;
import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.URL;

import javax.vecmath.Matrix4f;
import javax.vecmath.Point2f;
import javax.vecmath.Point3f;

import com.jhlabs.map.Ellipsoid;
import com.jhlabs.map.proj.MercatorProjection;
import com.jhlabs.map.proj.Projection;


public class GeneralProjProjection {

	protected int w;
	protected int h;

	private int minx, miny, maxx, maxy;

	protected int[][] vIndices;
	protected float[][] vertexWeights;

	private final Projection projection;

	public GeneralProjProjection(Projection proj) {
		this.projection = proj;
	}

	public static GeneralPath readDatFile(String file) throws IOException {
		return readDatFile(new BufferedReader(new FileReader(file)));
	}

	public static GeneralPath readDatFile(URL url) throws IOException {
	        BufferedReader in = new BufferedReader(
	        	new InputStreamReader(url.openStream()));
	        return readDatFile(in);
	}

	private static GeneralPath readDatFile(BufferedReader buf) throws IOException {
		GeneralPath path = new GeneralPath();
		String line;
		boolean closed = true;
		int lineno = 0;
		while((line = buf.readLine()) != null) {
			lineno++;
			line = line.trim();
			if(line.equals("")) continue;
			if(line.startsWith(">")) {
				closed = true;
				continue;
			}
			String[] toks = line.split("\\s");
			try {
				float x = Float.parseFloat(toks[0]);
				float y = Float.parseFloat(toks[1]);
				if(closed)
					path.moveTo(x, y);
				else
					path.lineTo(x, y);
			} catch(Exception e) {
				throw new RuntimeException("Error parsing line " + lineno + ": " + line, e);
			}
			closed = false;
		}
		buf.close();
		return path;
	}

	public int getWidth() {
		return w;
	}

	public int getHeight() {
		return h;
	}

	public static GeneralPath createLines() {
		GeneralPath lines = new GeneralPath();
		int pad = 0;
		int nLos = 24;
		int nLas = 500;
		for(int lo = 0; lo <= nLos; lo++) {
			boolean first = true;
			for(int la = pad; la <= nLas - pad; la++) {
				float x = (lo - nLos/2) * (360.0f / nLos);
				float y = (la - nLas/2) * (180.0f / nLas);
				if(first)
					lines.moveTo(x, y);
				else
					lines.lineTo(x, y);
				first = false;
			}
		}

		nLos = 500;
		nLas = 12;
		for(int la = 0; la <= nLas; la++) {
			boolean first = true;
			for(int lo = pad; lo <= nLos - pad; lo++) {
				float x = (lo - nLos/2) * (360.0f / nLos);
				float y = (la - nLas/2) * (180.0f / nLas);
				if(first)
					lines.moveTo(x, y);
				else
					lines.lineTo(x, y);
				first = false;
			}
		}
		return lines;
	}

	public GeneralPath createOutline() {
		GeneralPath outline = new GeneralPath();
		int nLos = 500;
		int nLas = 500;
		boolean first = true;

		for(int lo = 0; lo <= nLos; lo++) {
			float x = (lo - nLos / 2) * (360.0f / nLos);
			float y = -90;
			if(first)
				outline.moveTo(x, y);
			else
				outline.lineTo(x, y);
			first = false;
		}

		for(int la = 0; la <= nLos; la++) {
			float x = 180;
			float y = (la - nLas/2) * (180.0f / nLas);
			outline.lineTo(x, y);
		}

		for(int lo = nLos; lo >= 0; lo--) {
			float x = (lo - nLos / 2) * (360.0f / nLos);
			float y = 90;
			outline.lineTo(x, y);
		}

		for(int la = nLas; la >= 0; la--) {
			float x = -180;
			float y = (la - nLas/2) * (180.0f / nLas);
			outline.lineTo(x, y);
		}

		return outline;
	}

	private static double findRadiusForWidth(Projection projection, int w, double rLower, double rUpper) {

		double rGuess = (rUpper + rLower) / 2;

		if(rGuess == rLower || rGuess == rUpper)
			return rGuess;

		projection.setEllipsoid(new Ellipsoid("", rGuess, rGuess, 0.0D, ""));
		projection.initialize();

		int minx = Integer.MAX_VALUE;
		int maxx = -minx;
		int miny = minx;
		int maxy = -minx;

		Point2D.Double tmpin = new Point2D.Double();
		Point2D.Double tmpout = new Point2D.Double();
		double dlon = projection.getMaxLongitude() - projection.getMinLongitude();
		double dlat = projection.getMaxLatitude() - projection.getMinLatitude();

		for(int lo = 0; lo <= 100; lo++) {
			for(int la = 0; la <= 100; la++) {
				tmpin.x = (lo - 50) * dlon / 100;
				tmpin.y = (la - 50) * dlat / 100;
				tmpin.x = (lo - 50) * dlon / 100;
				tmpin.y = (la - 50) * dlat / 100;
				projection.transformRadians(tmpin, tmpout);
				if(tmpout.x < minx) minx = (int)tmpout.x;
				if(tmpout.y < miny) miny = (int)tmpout.y;
				if(tmpout.x > maxx) maxx = (int)tmpout.x;
				if(tmpout.y > maxy) maxy = (int)tmpout.y;
			}
		}
		int foundWidth = maxx - minx + 1;

		if(w == foundWidth)
			return rGuess;

		else if(foundWidth > w)
			return findRadiusForWidth(projection, w, rLower, rGuess);
		else
			return findRadiusForWidth(projection, w, rGuess, rUpper);
	}

	public void determineDimensions(int width) {
		double globeRadius = findRadiusForWidth(projection, width, 0, 1000f);
		projection.setEllipsoid(new Ellipsoid("", globeRadius, globeRadius, 0.0D, ""));
		projection.initialize();

		if(!projection.hasInverse())
			throw new RuntimeException("Inverse mapping missing");

		minx = Integer.MAX_VALUE;
		maxx = -minx;
		miny = minx;
		maxy = -minx;

		double dlon = projection.getMaxLongitude() - projection.getMinLongitude();
		double dlat = projection.getMaxLatitude() - projection.getMinLatitude();

		Point2D.Double tmpin = new Point2D.Double();
		Point2D.Double tmpout = new Point2D.Double();
		for(int lo = 0; lo <= 100; lo++) {
			for(int la = 0; la <= 100; la++) {
				tmpin.x = (lo - 50) * dlon / 100;
				tmpin.y = (la - 50) * dlat / 100;
				tmpin.x = (lo - 50) * dlon / 100;
				tmpin.y = (la - 50) * dlat / 100;
				projection.transformRadians(tmpin, tmpout);
				if(tmpout.x < minx) minx = (int)tmpout.x;
				if(tmpout.y < miny) miny = (int)tmpout.y;
				if(tmpout.x > maxx) maxx = (int)tmpout.x;
				if(tmpout.y > maxy) maxy = (int)tmpout.y;
			}
		}

		this.w = maxx - minx + 1;
		this.h = maxy - miny + 1;
	}

	public void prepareForProjection(SphericalMaxProjection smp, int width) {
		determineDimensions(width);
		vIndices = new int[w * h][3];
		vertexWeights = new float[w * h][3];
		Point3f p = new Point3f();
		Point3f[] vertices = smp.getSphere().getVertices();
		int index = 0;
		for(int y = maxy; y >= miny; y--) {
			for(int x = minx; x <= maxx; x++) {
				Point2D.Double in = new Point2D.Double(x, y);
				Point2D.Double out = new Point2D.Double();
				out = projection.inverseTransformRadians(in, out);

				double lon = out.x;
				double lat = out.y;

				if(lat >= projection.getMaxLatitude() ||
					lat <= projection.getMinLatitude() ||
					lon >= projection.getMaxLongitude() ||
					lon <= projection.getMinLongitude()) {
					index++;
					continue;
				}

				double sLon = (float)Math.sin(lon);
				double cLon = (float)Math.cos(lon);
				double sLat = (float)Math.sin(lat);
				double cLat = (float)Math.cos(lat);

				smp.getPoint(sLon, cLon, sLat, cLat, p);

				smp.getThreeNearestVertexIndices(p, vIndices[index]);

				// interpolate according to distance
				float d0 = 1 / p.distance(vertices[vIndices[index][0]]);
				float d1 = 1 / p.distance(vertices[vIndices[index][1]]);
				float d2 = 1 / p.distance(vertices[vIndices[index][2]]);
				float sum = d0 + d1 + d2;
				vertexWeights[index][0] = d0 / sum;
				vertexWeights[index][1] = d1 / sum;
				vertexWeights[index][2] = d2 / sum;

				index++;
			}
		}
	}

	public GeneralPath transform(GeneralPath in) {
		GeneralPath out = new GeneralPath();
		PathIterator it = in.getPathIterator(null);
		float[] seg = new float[6];

		Point2D.Double din = new Point2D.Double();
		Point2D.Double dout = new Point2D.Double();

		while(!it.isDone()) {
			int l = it.currentSegment(seg);
			din.x = seg[0];
			din.y = seg[1];
			if(din.x > projection.getMaxLongitudeDegrees() ||
					din.x < projection.getMinLongitudeDegrees() ||
					din.y > projection.getMaxLatitudeDegrees() ||
					din.y < projection.getMinLatitudeDegrees())
				l = PathIterator.SEG_MOVETO;
			projection.transform(din, dout);
			float x = (float)dout.x - minx;
			float y = maxy - (float)dout.y;

			if(!projection.inside(din.x, din.y))
				l = PathIterator.SEG_MOVETO;

			if(l == PathIterator.SEG_MOVETO) {
				out.moveTo(x, y);
			} else {
				out.lineTo(x, y);
			}
			it.next();
		}
		return out;
	}

	public void drawInto(ImageProcessor ip, int value, int linewidth, GeneralPath path) {
		ip.setValue(value);
		ip.setLineWidth(linewidth);
		PathIterator it = path.getPathIterator(null);
		float[] seg = new float[6];
		int px = 0, py = 0;
		while(!it.isDone()) {
			int l = it.currentSegment(seg);
			int x = Math.round(seg[0]);
			int y = Math.round(seg[1]);
			double d = Math.sqrt((px - x) * (px - x) + (py - y) * (py - y));

			if(l == PathIterator.SEG_MOVETO || d > 15)
				ip.moveTo(x, y);
			else
				ip.lineTo(x, y);
			px = x;
			py = y;
			it.next();
		}
	}

//	public static GeneralPath filter(GeneralPath in) {
//		GeneralPath out = new GeneralPath();
//		PathIterator it = in.getPathIterator(null);
//		float[] seg = new float[6];
//		float px = 0, py = 0;
//
//		int n = 0;
//		GeneralPath next = new GeneralPath();
//
//		while(!it.isDone()) {
//			int l = it.currentSegment(seg);
//
//			double d = Math.sqrt((px - seg[0]) * (px - seg[0]) + (py - seg[1]) * (py - seg[1]));
//			if(d > 20) {
//				l = PathIterator.SEG_MOVETO;
//			}
//
//			if(l == PathIterator.SEG_MOVETO) {
//				if(n > 50)
//					out.append(next.getPathIterator(null), false);
//				next = new GeneralPath();
//				next.moveTo(seg[0], seg[1]);
//				n = 0;
//			}
//			else {
//				next.lineTo(seg[0], seg[1]);
//				n++;
//			}
//			px = seg[0];
//			py = seg[1];
//			it.next();
//		}
//		if(n > 50)
//			out.append(next.getPathIterator(null), false);
//		return out;
//	}

	public int[] setSphereColors(SphericalMaxProjection smp, ColorProcessor ip) {
		Point2f polar = new Point2f();
		Point2D.Double out = new Point2D.Double();
		int[] colors = new int[smp.getSphere().nVertices];
		int idx = 0;
		for(Point3f v : smp.getSphere().getVertices()) {
			smp.getPolar(v, polar);
			projection.transformRadians(polar.x, polar.y, out);
			int x = (int)out.x - minx;
			int y = maxy - (int)out.y;

			colors[idx++] = ip.get(x, y);
		}
		return colors;
	}


	public static void savePath(GeneralPath path, String file, int w, int h, boolean fill) throws IOException {
		// path = filter(path);

		PrintWriter out = new PrintWriter(new FileWriter(file));
		out.println("%!PS-Adobe-3.0 EPSF-3.0");
		out.println("%%BoundingBox: 0 0 " + w + " " + h);
		out.println("%%EndComments");

		// Draw background:
//		out.println("0 0 moveto");
//		out.println(w + " 0 lineto");
//		out.println(w + " " + h + " lineto");
//		out.println("0 " + h + " lineto");
//		out.println("closepath");
//		out.println("0 setcolor");
//		out.println("fill");

		out.println("0.3 setlinewidth");

		PathIterator it = path.getPathIterator(null);
		float[] seg = new float[6];
		int px = 0, py = 0;
		while(!it.isDone()) {
			int l = it.currentSegment(seg);
			int x = Math.round(seg[0]);
			// int y = Math.round(seg[1]);
			int y = h - Math.round(seg[1]);
			double d = Math.sqrt((px - x) * (px - x) + (py - y) * (py - y));

			if(l == PathIterator.SEG_MOVETO
					|| d > 15
				)
				out.println(x + " " + y + " moveto");
			else
				out.println(x + " " + y + " lineto");
			px = x;
			py = y;
			it.next();
		}
		out.println("0 setcolor");
		if(fill)
			out.println("fill");
		else
			out.println("stroke");
		out.close();
	}

	// in needs to be in rad
	public GeneralPath transform(GeneralPath in, Matrix4f mat, SphericalMaxProjection smp) {
		GeneralPath out = new GeneralPath();
		PathIterator it = in.getPathIterator(null);
		float[] seg = new float[6];

		Point2f din = new Point2f();

		mat = new Matrix4f(mat);
		mat.invert();

		while(!it.isDone()) {
			int l = it.currentSegment(seg);
			din.x = (float)(seg[0] * Math.PI / 180.0);
			din.y = (float)(seg[1] * Math.PI / 180.0);

			smp.applyTransformToLonLat(din, mat);

			// transform back to degrees
			din.x = (float)(180 * din.x / Math.PI);
			din.y = (float)(180 * din.y / Math.PI);

			if(l == PathIterator.SEG_MOVETO) {
				out.moveTo(din.x, din.y);
			} else {
				out.lineTo(din.x, din.y);
			}
			it.next();
		}
		return out;
	}

	public ImageProcessor project(short[] maxima) {
		ShortProcessor ip = new ShortProcessor(w, h);
		for(int y = 0; y < h; y++) {
			for(int x = 0; x < w; x++) {
				int index = y * w + x;

				float v0 = vertexWeights[index][0] * (maxima[vIndices[index][0]] & 0xffff);
				float v1 = vertexWeights[index][1] * (maxima[vIndices[index][1]] & 0xffff);
				float v2 = vertexWeights[index][2] * (maxima[vIndices[index][2]] & 0xffff);

				ip.setf(x, y, (short)(v0 + v1 + v2));
			}
		}
		return ip;
	}

	public ImageProcessor projectColor(int[] maxima) {
		ColorProcessor ip = new ColorProcessor(w, h);
		for(int y = 0; y < h; y++) {
			for(int x = 0; x < w; x++) {
				int index = y * w + x;

				float w0 = vertexWeights[index][0];
				float w1 = vertexWeights[index][1];
				float w2 = vertexWeights[index][2];

				int m0 = maxima[vIndices[index][0]];
				int m1 = maxima[vIndices[index][1]];
				int m2 = maxima[vIndices[index][2]];

				int r = (int)(w0 * ((m0 & 0xff0000) >> 16) + w1 * ((m1 & 0xff0000) >> 16) + w2 * ((m2 & 0xff0000) >> 16));
				int g = (int)(w0 * ((m0 & 0xff00)   >>  8) + w1 * ((m1 & 0xff00)   >>  8) + w2 * ((m2 & 0xff00)   >>  8));
				int b = (int)(w0 * ((m0 & 0xff))           + w1 * ((m1 & 0xff))           + w2 * ((m2 & 0xff)));

				ip.set(x, y, (r << 16) + (g << 8) + b);
			}
		}
		return ip;
	}

	public static void main(String[] args) {
		GeneralProjProjection g = new GeneralProjProjection(new MercatorProjection());
		g.determineDimensions(576);
		// need h = 576
		System.out.println(g.w + " x " + g.h);
	}
}
