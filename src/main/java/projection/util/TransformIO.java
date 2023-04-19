package projection.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.scijava.vecmath.Matrix4f;

public class TransformIO {

	public static String toString(Matrix4f m) {
		if(m == null)
			return "null";
		StringBuilder sb = new StringBuilder();
		for(int r = 0; r < 4; r++)
			for(int c = 0; c < 4; c++)
				sb.append(m.getElement(r,  c)).append(" ");
		return sb.toString();
	}

	public static Matrix4f fromString(String s) {
		if(s.equals("null"))
			return null;
		String[] toks = s.split(" ");
		float[] matrix = new float[16];
		for(int i = 0; i < 16; i++)
			matrix[i] = Float.parseFloat(toks[i]);
		return new Matrix4f(matrix);
	}

	public static void saveTransformation(File file, Matrix4f m) throws IOException {
		PrintWriter out = new PrintWriter(new FileWriter(file));
		if(m == null)
			out.println("null");
		else {
			for(int r = 0; r < 4; r++) {
				for(int c = 0; c < 4; c++) {
					out.print(m.getElement(r,  c) + " ");
				}
			}
			out.println();
		}
		out.close();
	}

	public static Matrix4f loadTransformation(File file) throws IOException {
		BufferedReader in = new BufferedReader(new FileReader(file));
		String line = in.readLine();
		in.close();
		if(line.equals("null"))
			return null;
		String[] toks = line.split(" ");
		float[] matrix = new float[16];
		for(int i = 0; i < 16; i++)
			matrix[i] = Float.parseFloat(toks[i]);
		return new Matrix4f(matrix);
	}
}
