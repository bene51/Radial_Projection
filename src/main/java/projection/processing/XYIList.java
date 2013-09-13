package projection.processing;

public class XYIList {

	int[] x;
	int[] y;
	int[] i;

	int size = 0;

	public XYIList() {
		x = new int[10];
		y = new int[10];
		i = new int[10];
		size = 0;
	}

	public XYIList(int N) {
		x = new int[N];
		y = new int[N];
		i = new int[N];
		size = 0;
	}

	public int size() {
		return size;
	}

	private void ensureSize(int n) {
		if(x.length >= n)
			return;
		int[] tmp = new int[n];
		System.arraycopy(x, 0, tmp, 0, x.length);
		x = tmp;
		tmp = new int[n];
		System.arraycopy(y, 0, tmp, 0, y.length);
		y = tmp;
		tmp = new int[n];
		System.arraycopy(i, 0, tmp, 0, i.length);
		i = tmp;
	}

	public final void add(int xn, int yn, int in) {
		if(x.length == size)
			ensureSize(x.length + x.length / 2);
		x[size] = xn;
		y[size] = yn;
		i[size] = in;
		size++;
	}

	public static final XYIList addAll(XYIList[] lists) {
		int L = 0;
		for(XYIList list : lists)
			L += list.size();
		XYIList newl = new XYIList(L);
		int offs = 0;
		for(XYIList list : lists) {
			int n = list.size();
			System.arraycopy(list.x, 0, newl.x, offs, n);
			System.arraycopy(list.y, 0, newl.y, offs, n);
			System.arraycopy(list.i, 0, newl.i, offs, n);
			offs += n;
		}
		newl.size = L;
		return newl;
	}
}
