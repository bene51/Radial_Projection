package projection.processing;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.vecmath.Point3f;
import javax.vecmath.Vector3f;

public class MaximumProjector {
	// These fields are set in prepareForProjection();
	private int[][] lutxy;
	private int[][] luti;
	private int[][] luts;

	private short[][] maxima;

	private int[] usedVertexIndices;
	private byte[] forSavingVertices;

	private final SphericalMaxProjection smp;

	public MaximumProjector(SphericalMaxProjection smp) {
		this.smp = smp;
	}

	/*
	 * z starts with 0;
	 */
	public void projectPlaneMultilayer(int z, short[] ip) {
		for(int i = 0; i < luti[z].length; i++) {
			float v = ip[lutxy[z][i]] & 0xffff;
			if(v > (maxima[luts[z][i]][luti[z][i]] & 0xffff))
				maxima[luts[z][i]][luti[z][i]] = (short)v;
		}
	}

	public void resetMaxima() {
		for(int s = 0; s < maxima.length; s++)
			for(int v = 0; v < maxima[s].length; v++)
				maxima[s][v] = 0;
	}

	public void saveVertices(int l, File file) throws IOException {
		int i = 0;
		for(short pixel : maxima[l]) {
			forSavingVertices[i++] = (byte) pixel;
			forSavingVertices[i++] = (byte) (pixel >> 8);
		}
		File tmp = new File(file.getAbsolutePath() + ".part");
		FileOutputStream fos = new FileOutputStream(tmp);
		fos.write(forSavingVertices);
		fos.close();
		tmp.renameTo(file);
	}

	public static void loadVertices(File file, int[] usedVertexIndices, short[] maxima) throws IOException {
		byte[] buffer = new byte[2 * usedVertexIndices.length];
		FileInputStream in = new FileInputStream(file);
		in.read(buffer);
		in.close();
		int i = 0;
		for(int vIdx : usedVertexIndices) {
			int low  = 0xff & buffer[i++];
			int high = 0xff & buffer[i++];
			maxima[vIdx] = (short)((high << 8) | low);
		}
	}

	public void saveVertexIndices(File file) throws IOException {
		PrintStream out = new PrintStream(file);
		for(int v : usedVertexIndices)
			out.println(Integer.toString(v));
		out.close();
	}

	public static int[] loadVertexIndices(File file) throws IOException {
		ArrayList<Integer> vs = new ArrayList<Integer>();
		BufferedReader in = new BufferedReader(new FileReader(file));
		String line;
		while((line = in.readLine()) != null)
			vs.add(Integer.parseInt(line));
		int[] ret = new int[vs.size()];
		for(int i = 0; i < ret.length; i++)
			ret[i] = vs.get(i);
		return ret;
	}

	public void prepareForProjection(final int w, final int h, final int d, final double pw, final double ph, final double pd, final double layerWidth, final int nLayers, final FusionWeight weighter) {

		long totalStart = System.currentTimeMillis();
		final int nProcessors = Runtime.getRuntime().availableProcessors();
		ExecutorService exec = Executors.newFixedThreadPool(nProcessors);

		final XYIList[][] correspondences = new XYIList[d][nProcessors];
		final XYIList[] all_correspondences = new XYIList[d];

		final double t2 = nLayers * layerWidth / 2;

		final ArrayList<Integer> allValidVertexIndices = new ArrayList<Integer>();
		final IndexedTriangleMesh sphere = smp.getSphere();
		final Point3f center = smp.getCenter();
		for(int proc = 0; proc < nProcessors; proc++) {
			final int currentProc = proc;
			exec.execute(new Runnable() {
				@Override
				public void run() {
					int nVerticesPerThread = (int)Math.ceil(sphere.nVertices / (double)nProcessors);
					int startV = currentProc * nVerticesPerThread;
					int lenV = Math.min((currentProc + 1) * nVerticesPerThread, sphere.nVertices);
					ArrayList<Integer> validVertexIndices = new ArrayList<Integer>();
					for(int vIndex = startV; vIndex < lenV; vIndex++) {
						Point3f vertex = sphere.getVertices()[vIndex];
						float weight = weighter.getWeight(vertex.x, vertex.y, vertex.z);
						if(weight != 0)
							validVertexIndices.add(vIndex);
					}
					synchronized(MaximumProjector.this) {
						allValidVertexIndices.addAll(validVertexIndices);
					}
				}
			});
		}
		try {
			exec.shutdown();
			exec.awaitTermination(30, TimeUnit.MINUTES);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		usedVertexIndices = new int[allValidVertexIndices.size()];
		forSavingVertices = new byte[2 * usedVertexIndices.length];
		for(int i = 0; i < usedVertexIndices.length; i++)
			usedVertexIndices[i] = allValidVertexIndices.get(i);

		exec = Executors.newFixedThreadPool(nProcessors);
		for(int proc = 0; proc < nProcessors; proc++) {
			final int currentProc = proc;
			exec.execute(new Runnable() {
				@Override
				public void run() {
					Vector3f dx = new Vector3f();
					Point3f pos = new Point3f();

					for(int i = 0; i < d; i++)
						correspondences[i][currentProc] = new XYIList();

					int nVerticesPerThread = (int)Math.ceil(allValidVertexIndices.size() / (double)nProcessors);
					int startV = currentProc * nVerticesPerThread;
					int lenV = Math.min((currentProc + 1) * nVerticesPerThread, allValidVertexIndices.size());
					for(int vvi = startV; vvi < lenV; vvi++) {
						int vIndex = allValidVertexIndices.get(vvi);
						Point3f vertex = sphere.getVertices()[vIndex];
						float weight = weighter.getWeight(vertex.x, vertex.y, vertex.z);
						if(weight == 0) {
							System.out.println("weight is 0. Shouldn't");
							continue;
						}

						dx.sub(vertex, center);
						dx.normalize();

						// calculate the distance needed to move to the neighbor pixel
						double dmax = Math.max(Math.abs(dx.x / pw), Math.max(Math.abs(dx.y / ph), Math.abs(dx.z / pd)));
						double ddx = dx.x / dmax;
						double ddy = dx.y / dmax;
						double ddz = dx.z / dmax;
						double scale = Math.sqrt(ddx * ddx + ddy * ddy + ddz * ddz);

						for(int s = 0; s < nLayers; s++) {
							double start = s * layerWidth - t2;
							double end = start + layerWidth;

							for(double i = start; i <= end; i += scale) {
								pos.scaleAdd((float)i, dx, vertex);

								// nearest neighbor in pixel dims
								int x = (int)Math.round(pos.x / pw);
								int y = (int)Math.round(pos.y / ph);
								int z = (int)Math.round(pos.z / pd);

								// only add it if the pixel is inside the image
								if(x >= 0 && x < w && y >= 0 && y < h && z >= 0 && z < d)
									correspondences[z][currentProc].add(y * w + x, s, vvi);
							}
						}
					}
					System.out.println("thread " + currentProc + " done");
				}
			});
		}

		try {
			exec.shutdown();
			exec.awaitTermination(30, TimeUnit.MINUTES);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}

		lutxy = new int[d][];
		luti = new int[d][];
		luts = new int[d][];
		for(int i = 0; i < d; i++) {
			all_correspondences[i] = XYIList.addAll(correspondences[i]);
			lutxy[i] = all_correspondences[i].x;
			luts[i]  = all_correspondences[i].y;
			luti[i]  = all_correspondences[i].i;
		}

		long totalEnd = System.currentTimeMillis();

		maxima = new short[nLayers][usedVertexIndices.length];
		System.out.println("Overall time " + (totalEnd - totalStart) + " ms");
	}
}
