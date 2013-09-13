package projection.processing;

import ij.IJ;
import ij.ImageJ;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import fiji.util.gui.GenericDialogPlus;


public class Tmp {

	private int n = 3;
	private GenericDialogPlus gd;

	void showDialog() {
		gd = new GenericDialogPlus("Add directory");

		gd.addButton("Add directory", new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				IJ.log("button pressed");
				gd.dispose();
				n++;
				showDialog();
			}
		});
		for(int i = 0; i < n; i++)
			gd.addDirectoryField("dir " + i, "");

		gd.showDialog();
	}

	public static void main(String[] args) {
		new ImageJ();
		new Tmp().showDialog();
	}
}
