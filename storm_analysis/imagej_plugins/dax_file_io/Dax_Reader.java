//ImageJ reader plugin for Zhuang lab .dax files
//Translated from the Python versions in storm-analysis/sa_library and 
//based on sample IJ reader/writer plugins by Albert Cardona at 
//http://albert.rierol.net/imagej_programming_tutorials.html.
//Evan Heller, 10/2015

import java.io.*;  
import ij.*;  
import ij.io.*;
import ij.plugin.PlugIn;
import ij.plugin.FileInfoVirtualStack;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

public class Dax_Reader extends ImagePlus implements PlugIn {

    public void run(String arg) {  
        String path = getPath(arg);  
        if (null == path) return;  
        if (!parse(path)) return;  
        if (null == arg || 0 == arg.trim().length()) this.show(); // was opened by direct call to the plugin  
        // not via HandleExtraFileTypes which would  
        // have given a non-null arg.  
    }  

    private String getPath(String arg) {  
        if (null != arg) {  
            if (0 == arg.indexOf("http://")  || new File(arg).exists()) return arg;  
        }  
        // else, ask:  
        OpenDialog od = new OpenDialog("Choose a .dax file", null);  
        String dir = od.getDirectory();  
        if (null == dir) return null; // dialog was canceled  
        dir = dir.replace('\\', '/'); // Windows safe  
        if (!dir.endsWith("/")) dir += "/";  
        return dir + od.getFileName();  
    }  

    private boolean parse(String path) {  
        File mydax = new File(path);
        String dirname = mydax.getParent();
        String filename = mydax.getName();
        if (dirname.length() > 0) dirname += "/";

        String inf_file = dirname + filename.split("\\.")[0] + ".inf";

        //Image parameters
        int image_height = 256;
        int image_width = 256;
        int number_frames = 0, scalemax=0, scalemin=0;
        boolean bigendian = false;
        float lock_target, stage_x = 0, stage_y = 0;

        //Extract the movie information from the associated inf file
        Pattern size_re = Pattern.compile("frame dimensions = ([\\d]+) x ([\\d]+)");
        Pattern length_re = Pattern.compile("number of frames = ([\\d]+)");
        Pattern endian_re = Pattern.compile(" (big|little) endian");
        Pattern stagex_re = Pattern.compile("Stage X = ([\\d\\.\\-]+)");
        Pattern stagey_re = Pattern.compile("Stage Y = ([\\d\\.\\-]+)");
        Pattern lock_target_re = Pattern.compile("Lock Target = ([\\d\\.\\-]+)");
        Pattern scalemax_re = Pattern.compile("scalemax = ([\\d\\.\\-]+)");
        Pattern scalemin_re = Pattern.compile("scalemin = ([\\d\\.\\-]+)");

        try {
            BufferedReader br = new BufferedReader(new FileReader(inf_file));
            String line;
            Matcher m;

            while ( (line = br.readLine()) != null) {
                m = size_re.matcher(line);
                if (m.find()) {
                    image_width = Integer.parseInt(m.group(1));
                    image_height = Integer.parseInt(m.group(2));
                }

                m = length_re.matcher(line);
                if ( m.find() ) number_frames = Integer.parseInt(m.group(1));

                m = endian_re.matcher(line);
                if ( m.find() ) { 
                    if (m.group(1).equals("big")) bigendian = true;
                }

                m = stagex_re.matcher(line);
                if ( m.find() ) stage_x = Float.parseFloat(m.group(1));

                m = stagey_re.matcher(line);
                if ( m.find() ) stage_y = Float.parseFloat(m.group(1));

                m = lock_target_re.matcher(line);
                if (m.find() ) lock_target = Float.parseFloat(m.group(1));

                m = scalemax_re.matcher(line);
                if (m.find() ) scalemax = Integer.parseInt(m.group(1));

                m = scalemin_re.matcher(line);
                if ( m.find() ) scalemin = Integer.parseInt(m.group(1));
            }

            br.close();

        } 
        catch (IOException e) {
            System.err.println("Caught IOException: " + e.getMessage());
            return false;
        }

        FileInfo fi = new FileInfo();
        fi.fileType=2;
        fi.fileFormat=fi.TIFF;
        fi.directory=dirname;
        fi.fileName=filename.split("\\.")[0] + ".dax";
        fi.width=image_width;
        fi.height=image_height;
        fi.nImages=number_frames;
        fi.gapBetweenImages = 0;
        fi.intelByteOrder = !bigendian;
        fi.whiteIsZero = false;
        fi.longOffset = fi.offset = 0;

        try {
            ImagePlus imp;

            //Open regular stack if small enough; otherwise, virtual
            if (number_frames <= 500) {
                FileOpener fo = new FileOpener(fi);  
                imp = fo.open(false);  
            } 
            else {
                FileInfoVirtualStack fv = new FileInfoVirtualStack(fi, false);
                imp = new ImagePlus("",fv);
            }

            this.setStack(imp.getTitle(), imp.getStack());  
            this.setCalibration(imp.getCalibration());  
            Object obinfo = imp.getProperty("Info");  
            if (null != obinfo) this.setProperty("Info", obinfo);  
            this.setFileInfo(imp.getOriginalFileInfo());  

        }
        catch (Exception e) {
            e.printStackTrace(); 
            return false;
        }

        return true;
    }
}
