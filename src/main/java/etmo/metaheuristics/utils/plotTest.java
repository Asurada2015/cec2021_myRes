package etmo.metaheuristics.utils;

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

public class plotTest {

    double[][] pf1 , pf2;
    MetricsUtil util ;

    public plotTest() {
        util = new MetricsUtil();
        pf1 = util.readFront("PF/cec2019/concave.pf");
        pf2 = util.readFront("PF/cec2019/convex.pf");
    }

    public void doTest(String name, String title){
        XYSeriesCollection xyDataset = new XYSeriesCollection();
        XYSeries series1 = new XYSeries("pf1");
        for (int i = 0; i < pf1.length; i++){
            series1.add(pf1[i][0],pf1[i][1]);
        }
        xyDataset.addSeries(series1);

        XYSeries series2 = new XYSeries("pf2");
        for (int i = 0; i < pf2.length; i++){
            series2.add(pf2[i][0],pf2[i][1]);
        }
        xyDataset.addSeries(series2);

        JFreeChart chart = ChartFactory.createScatterPlot("test","x","y",xyDataset);

        ChartFrame frame = new ChartFrame("2D scatter plot", chart, true);
        XYPlot xyplot = (XYPlot) chart.getPlot();
        xyplot.setBackgroundPaint(Color.white);//设置背景面板颜色
        ValueAxis vaaxis = xyplot.getDomainAxis();
        vaaxis.setAxisLineStroke(new BasicStroke(1.5f));//设置坐标轴粗细

        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.pack();
        frame.setVisible(true);

    }

    public static void main(String[] args) {
        plotTest test2d = new plotTest();
        test2d.doTest("test2d", "2objective");
    }
}
