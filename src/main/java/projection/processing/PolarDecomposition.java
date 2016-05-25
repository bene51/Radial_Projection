package projection.processing;

import ij.gui.Plot;

import java.awt.Color;

import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Matrix4d;
import javax.vecmath.Vector3d;

public class PolarDecomposition {

	public final double rotationAxisX, rotationAxisY, rotationAxisZ;
	public final double angle; // in degree
	public final double scaleX, scaleY, scaleZ;
	public final double shiftX, shiftY, shiftZ;
	public final double determinant;

	public PolarDecomposition(Matrix4d m) {
		Vector3d trans = new Vector3d();
		m.get(trans);
		shiftX = trans.x;
		shiftY = trans.y;
		shiftZ = trans.z;

		Matrix3d q = new Matrix3d();
		Matrix3d aff = new Matrix3d();
		m.getRotationScale(aff);
		determinant = polarDecomposition(aff, q);

		AxisAngle4d aa = new AxisAngle4d();
		aa.set(q);
		rotationAxisX = aa.x;
		rotationAxisY = aa.y;
		rotationAxisZ = aa.z;
		angle = aa.getAngle() * 180 / Math.PI;

		Matrix4d ss = new Matrix4d(m);
		ss.setRotationScale(q);
		ss.invert();
		ss.mul(m);

		scaleX = ss.getM00();
		scaleY = ss.getM11();
		scaleZ = ss.getM22();
	}

	public static double frobenius(Matrix3d m1, Matrix3d m2) {
		double sum = 0;
		for(int i = 0; i < 3; i++) {
			for(int j = 0; j < 3; j++) {
				double d = m1.getElement(i, j) - m2.getElement(i, j);
				sum += d * d;
			}
		}
		return sum;
	}

	public static double polarDecomposition(Matrix3d m, Matrix3d q) {
		q.set(m);
		Matrix3d next = new Matrix3d();
		for(int i = 0; i < 100; i++) {
			next.set(q);
			next.transpose();
			next.invert();
			next.add(q);
			next.mul(0.5);
			if(frobenius(q, next) < 1e-10) {
				q.set(next);
				break;
			}
			q.set(next);
		}
		return q.determinant();
	}

	public static void main(String[] args) {
		double[] x = new double[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
		double[] y = new double[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
		Plot p = new Plot("bla", "x", "y", x, y);
		p.draw();
		p.setColor(Color.RED);
		double[] z = new double[] {5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
		p.addPoints(x, z, Plot.LINE);
		p.show();
	}
}
