package etmo.util.comparators;

import etmo.core.Solution;

import java.util.Comparator;

public class PdfComparator implements Comparator {
    @Override
    public int compare(Object o1, Object o2) {
        double pdf1 = ((Solution)o1).getPdf();
        double pdf2 = ((Solution)o2).getPdf();
        if (pdf1 > pdf2){
            return -1;
        }
        if (pdf1 < pdf2){
            return 1;
        }
        return 0;
    }
}
