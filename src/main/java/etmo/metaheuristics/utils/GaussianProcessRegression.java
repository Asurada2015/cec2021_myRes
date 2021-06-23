package etmo.metaheuristics.utils;

import smile.math.kernel.MercerKernel;

public class GaussianProcessRegression extends smile.regression.GaussianProcessRegression {
    public GaussianProcessRegression(MercerKernel kernel, Object[] regressors, double[] weight, double noise) {
        super(kernel, regressors, weight, noise);
    }


}
