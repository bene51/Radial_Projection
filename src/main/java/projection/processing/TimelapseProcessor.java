package projection.processing;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

public abstract class TimelapseProcessor {

	public abstract void processTimepoint(int tp) throws Exception;

	public void processParallel(int nTimepoints) {
		final int nProcessors = Runtime.getRuntime().availableProcessors();
		ExecutorService exec = Executors.newFixedThreadPool(nProcessors);
		int nFilesPerThread = (int)Math.ceil(nTimepoints / (double)nProcessors);

		for(int p = 0; p < nProcessors; p++) {
			final int start = p * nFilesPerThread;
			final int end = Math.min(nTimepoints, (p + 1) * nFilesPerThread);

			exec.submit(new Runnable() {
				@Override
				public void run() {
					for(int tp = start; tp < end; tp++) {
						try {
							processTimepoint(tp);
						} catch(Exception e) {
							System.err.println("Cannot process timepoint " + tp);
							e.printStackTrace();
						}
					}
				}
			});
		}
		try {
			exec.shutdown();
			exec.awaitTermination(300, TimeUnit.MINUTES);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}

	public void processSequential(int nTimepoints) {
		for(int tp = 0; tp < nTimepoints; tp++) {
			try {
				processTimepoint(tp);
			} catch(Exception e) {
				System.err.println("Cannot process timepoint " + tp);
				e.printStackTrace();
			}
		}
	}
}
