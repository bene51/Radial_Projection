package projection.processing;

import pal.math.ConjugateDirectionSearch;
import pal.math.MultivariateFunction;

public abstract class Optimizer {

	/**
	 * Real-world (unnormalized) parameters.
	 * @param parameters
	 * @return
	 */
	public abstract double calculateDifference(double[] parameters);

	public double optimize(double tol, double[] initial, double[] max) {
		ConjugateDirectionSearch CG = new ConjugateDirectionSearch();

		Refinement refinement = new Refinement(initial, max);
		double[] parameters = new double[refinement.getNumArguments()];
		double[] optimized;
		double badness = Float.MAX_VALUE;

		double maxInt = max(max);

		do {
			CG.optimize(refinement, parameters, tol,  tol);
			parameters = refinement.best;
			badness = refinement.evaluate(parameters);
			optimized = refinement.getParameters(parameters);
		} while(refinement.adjustInitial(parameters) > maxInt / 8);

		System.out.println("Difference = " + badness);
		System.arraycopy(optimized, 0, initial, 0, optimized.length);
		return refinement.bestValue;
	}

	private static double max(double[] array) {
		double max = array[0];
		for(int i = 1; i < array.length; i++)
			if(array[i] > max)
				max = array[i];
		return max;
	}

	private class Refinement implements MultivariateFunction {

		private final int N;

		private final double[] max;
		private final double[] factor;
		private double[] initial;

		private double bestValue = Double.MAX_VALUE;
		private double[] best;

		public Refinement(double[] initial, double[] max) {
			this.N = initial.length;
			this.initial = initial;
			this.max = max;
			this.factor = new double[N];
			double tmp = max(max);
			for(int i = 0; i < N; i++)
				factor[i] = max[i] / tmp;
		}

		/*
		 * @implements getNumArguments() in MultivariateFunction
		 */
		@Override
		public int getNumArguments() {
			return N;
		}

		/*
		 * @implements getLowerBound() in MultivariateFunction
		 */
		@Override
		public double getLowerBound(int n) {
			return -max[n] / factor[n];
		}

		/*
		 * @implements getUpperBound() in MultivariateFunction
		 */
		@Override
		public double getUpperBound(int n) {
			return max[n] / factor[n];
		}

		/*
		 * @implements evaluate() in MultivariateFunction
		 */
		@Override
		public double evaluate(double[] a) {
			double diff = calculateDifference(getParameters(a));
			if(diff < bestValue) {
				bestValue = diff;
				best = a.clone();
			}
			return diff;
		}

		/*
		 * Adjusts the initial parameters, resets the given array and
		 * returns the maximum absolute (normalized) adjustment.
		 * TODO why the normalized adjustment?
		 */
		public double adjustInitial(double[] a) {
			for(int i = 0; i < N; i++)
				initial[i] += a[i] * factor[i];

			double maxAdjust = 0;
			for (int i = 0; i < N; i++) {
				if (Math.abs(a[i]) > maxAdjust)
					maxAdjust = Math.abs(a[i]);
				a[i] = 0;
			}
			return maxAdjust;
		}

		/*
		 * Takes the optimized (normalized) 6-el parameter array
		 * and calculates the 9 real-world euler parameters.
		 */
		public double[] getParameters(double[] a) {
			double[] param = new double[N];
			for(int i = 0; i < N; i++)
				param[i] = a[i] * factor[i] + initial[i];
			return param;
		}
	}
}
