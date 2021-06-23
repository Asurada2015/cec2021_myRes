package etmo.metaheuristics.utils;


import smile.math.kernel.GaussianKernel;
import smile.regression.GaussianProcessRegression;

public class TestGP {
    public static void main(String[] args) {
        double[][] x = new double[10][2];
        double[] y = new double[10];
        for (int i = 1; i <= 10; i++){
            x[i - 1][0] = i;
            x[i - 1][1] = i + 1;
            y[i - 1] = x[i - 1][0] + x[i - 1][1];
        }
        GaussianProcessRegression<double[]> gpModel = GaussianProcessRegression.<double[]>fit(x, y, new GaussianKernel(0.6), 0.00);
//        System.out.println(gpModel);
//        System.out.println(gpModel.predict(new double[]{2.0,3.0}));

        String a = "22" + 3;
        System.out.println(a);

    }
}
