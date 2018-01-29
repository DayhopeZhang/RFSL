// Copyright (c) Ronghua Liang and Xiaorui Jiang
// Scientiﬁc ranking over heterogeneous academic hypernetwork. In AAAI, pages 20–26, 2016.

package utils;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class Utility {
    public static boolean isZero(double val) {
        return Math.abs(val) < Constants.EPSILON;
    }

    public static boolean isNumeric(String str) {
        Pattern pattern = Pattern.compile("[0-9]*");
        Matcher isNum = pattern.matcher(str);
        if (!isNum.matches()) {
            return false;
        }
        return true;
    }

    public static String ZeroPad( int number, int width ) {
        StringBuffer result = new StringBuffer("");
        for( int i = 0; i < width-Integer.toString(number).length(); i++ )
            result.append( "0" );
        result.append( Integer.toString(number) );
        return result.toString();
    }

    public static String ZeroPad( String number, int width ) {
        number = number.trim();
        StringBuffer result = new StringBuffer("");
        for( int i = 0; i < width-number.length(); i++ )
            result.append( "0" );
        result.append(number);
        return result.toString();
    }

    public static String padOSPath(String path) {
        String res = path;
        int sepLen = Constants.sep.length();
        int idx = res.lastIndexOf(Constants.sep);
        if (idx != res.length() - sepLen)
            res += Constants.sep;
        return res;
    }

    public static byte[] readFile(String filePath) throws IOException {

        File file =new File(filePath);
        if(filePath==null || filePath.equals(""))
        {
            throw new NullPointerException("Invalid File Path!");
        }
        long len = file.length();
        byte[] bytes = new byte[(int)len];

        BufferedInputStream bufferedInputStream=new BufferedInputStream(new FileInputStream(file));
        int r = bufferedInputStream.read( bytes );
        if (r != len)
            throw new IOException("Error Reading File!");
        bufferedInputStream.close();

        return bytes;

    }

    public static void writeFile(byte[] data, String filePath) throws IOException {
        File file =new File(filePath);
        file.getParentFile().mkdirs();
        BufferedOutputStream bufferedOutputStream=new BufferedOutputStream(new FileOutputStream(file));
        bufferedOutputStream.write(data);
        bufferedOutputStream.close();

    }

    public static void copyFile(String fromPath, String toPath) throws IOException {
        byte[] bytes = readFile(fromPath);
        writeFile(bytes, toPath);
    }

    public static void deleteFile(String path) {
        File f = new File(path);
        if (f.exists())
            f.delete();
    }

    /**
     *
     * @param real
     * @param precision
     * @return
     */
    public static String round(double real, int precision) {
        String formatString = "0.";
        for (int i = 0; i < precision; ++i)
            formatString += "0";
        DecimalFormat df = new DecimalFormat(formatString);		// keeping the precision
        return df.format(real);
//		df.setMaximumFractionDigits(precision);
//        // RoundingMode.DOWN v.s. RoundingMode.DOWN
//        nf.setRoundingMode(RoundingMode.UP);
//        return ZeroPad(nf.format(real), precision);
    }
}
