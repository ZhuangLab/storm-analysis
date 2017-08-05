//ImageJ writer plugin for Zhuang lab .dax files
//Translated from the Python versions in storm-analysis/sa_library and 
//based on sample IJ reader/writer plugins by Albert Cardona at 
//http://albert.rierol.net/imagej_programming_tutorials.html.
//Evan Heller, 10/2015

import java.io.*;  
import ij.*;  
import ij.io.*;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import ij.process.ImageConverter;

public class Dax_Writer implements PlugIn {  

    public void run(String arg) {  
        ImagePlus imp = WindowManager.getCurrentImage();  
        if (null == imp) return;  
        SaveDialog sd = new SaveDialog("Save Dax", "untitled", null);  
        String dir = sd.getDirectory();  
        if (null == dir) return; // user canceled dialog  
        dir = dir.replace('\\', '/'); // Windows safe  
        if (!dir.endsWith("/")) dir += "/";  
        saveDax(imp, dir + sd.getFileName());  
    }  

    static public void saveDax(ImagePlus imp, String path) {  
        File file = new File(path);  
        DataOutputStream dos = null;  

        try {  
            dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(file)));  

            // read data:  
            FileInfo fi = imp.getFileInfo(); 
            String inf_file = path.split("\\.")[0] + ".inf";

            FileWriter fw = new FileWriter(inf_file); 

            //Write a rudimentary .inf file 
            fw.write("binning = 1 x 1\n");
            fw.write("data type = 16 bit integers (binary, big endian)\n");
            fw.write("frame dimensions = " + fi.width + " x " + fi.height + "\n");
            fw.write("number of frames = " + fi.nImages + "\n");
            fw.write("Lock Target = 0.0\n");
            fw.write("x_start = 1\n");
            fw.write("x_end = " + fi.width + "\n");
            fw.write("y_start = 1\n");
            fw.write("y_end = " + fi.height + "\n");
            fw.flush();
            fw.close();

            ImageProcessor ip;

            //Covert to 16-bit grayscale for writing to .dax
            if (imp.getBitDepth() != 16) {
                ImageConverter ic = new ImageConverter(imp);
                ic.convertToGray16();
            }

            //Write the data to file
            for (int i=0; i<imp.getNSlices(); i++) {  
                imp.setSlice(i);
                ip= imp.getProcessor(); //Get ip of current slice
                short[] px = (short[])ip.getPixels();

                for (int j=0; j<px.length; j++) 
                    dos.writeShort(px[j]);

            }  
            dos.flush();  
            dos.close(); 

        } catch (Exception e) {  
            e.printStackTrace();  
        } 
    }  
}  
