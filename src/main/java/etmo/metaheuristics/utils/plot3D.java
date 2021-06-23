package etmo.metaheuristics.utils;

import tanling.matplot_4j.d3d.base.speg.Point3D;
import tanling.matplot_4j.d3d.facade.MatPlot3DMgr;

import java.util.ArrayList;

public class plot3D {
    public static void main(String[] args) {
        MatPlot3DMgr mgr = new MatPlot3DMgr();

        mgr.setDataInputType(MatPlot3DMgr.DATA_TYPE_DOTS);

        //*************************************************************//
        //Add your data here
        Point3D a = new Point3D(10, 1, 1);
        Point3D b = new Point3D(2, 20, 2);
        Point3D c=new Point3D(3,3,30);

        ArrayList<Point3D> aa = new ArrayList<>();
        aa.add(a);
        aa.add(b);
        aa.add(c);

        mgr.addData("Item 1", aa);




        mgr.setScaleX(1.2);
        mgr.setScaleY(1.2);
        mgr.setScaleZ(1.2);

        mgr.setSeeta(0.6);
        mgr.setBeita(1.0);

        mgr.setTitle("Demo");

        mgr.getProcessor().setCoordinateSysShowType(mgr.getProcessor().COORDINATE_SYS_ALWAYS_FURTHER);

        mgr.show();

    }
}
