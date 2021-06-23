package etmo.metaheuristics.utils;

import etmo.core.SolutionSet;
import etmo.qualityIndicator.util.MetricsUtil;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javax.swing.*;
import java.awt.*;

public class Plot2D_demo1 {
    double[][] pf;
    MetricsUtil util;

    public Plot2D_demo1(String pfName) {
        this.util = new MetricsUtil();
        this.pf = util.readFront(pfName);
    }

    public void doPlot2D(SolutionSet pop) throws InterruptedException {
        XYSeriesCollection xyDataset = new XYSeriesCollection();
        XYSeries pfSeries = new XYSeries("pf");
        for (int i = 0; i < pf.length; i++){
            pfSeries.add(pf[i][0],pf[i][1]);
        }
        xyDataset.addSeries(pfSeries);

        XYSeries popSeries = new XYSeries("solutions");
        for (int i = 0; i < pop.size(); i++){
            popSeries.add(pop.get(i).getObjective(0), pop.get(i).getObjective(1));
        }
        xyDataset.addSeries(popSeries);



        JFreeChart chart = ChartFactory.createScatterPlot("test","x","y",xyDataset);

        ChartFrame frame = new ChartFrame("2D scatter plot", chart, true);
        XYPlot xyplot = (XYPlot) chart.getPlot();


        xyplot.setBackgroundPaint(Color.white);//设置背景面板颜色
        ValueAxis vaaxis = xyplot.getDomainAxis();
        vaaxis.setAxisLineStroke(new BasicStroke(1.5f));//设置坐标轴粗细

        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setVisible(true);



//        Thread.sleep(1000);
//
//        xyDataset.removeSeries(pfSeries);
//        double[][] pf2 = util.readFront("PF/cec2019/concave.pf");
//        XYSeries pfSeries2 = new XYSeries("pf2");
//        for (int i = 0; i < pf2.length; i++){
//            pfSeries2.add(pf2[i][0],pf2[i][1]);
//        }
//        xyDataset.addSeries(pfSeries2);






    }
}
