package etmo.metaheuristics.utils;

import Jama.Matrix;
import etmo.qualityIndicator.util.MetricsUtil;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.DefaultXYDataset;

import javax.swing.*;
import java.awt.*;

public class plot2D {

    MetricsUtil util ;
    double[][] pf;
    double[][] transPf ;

    public plot2D() {
        util = new MetricsUtil();
        pf = util.readFront("PF/cec2019/concave.pf");
        transPf =  new double[2][pf.length];
        for (int i = 0; i < pf.length; i++){
            transPf[0][i] = pf[i][0];
            transPf[1][i] = pf[i][1];
        }
    }


    void doPlot_2d(double[][] data, String name, String title){
        DefaultXYDataset xydataset = new DefaultXYDataset();
        xydataset.addSeries(title, data);//设置点的图标title一般表示这画的是决策变量还是目标函数值
        JFreeChart chart = ChartFactory.createScatterPlot(name, "X", "Y", xydataset,
                PlotOrientation.VERTICAL, true, true, false);//设置表头，x轴，y轴，name表示问题的类型

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
        plot2D test2d = new plot2D();
        test2d.doPlot_2d(test2d.transPf,  "test2d", "2objective");
    }
}
