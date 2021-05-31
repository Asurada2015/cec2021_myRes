package etmo.metaheuristics.autoencoder;

import etmo.core.Solution;
import etmo.core.SolutionSet;
import etmo.util.JMException;
import org.deeplearning4j.nn.conf.layers.AutoEncoder;

import java.util.Arrays;

public class Utils {
    public static double[][] MappingViaAE(SolutionSet pop, int task1, int task2, int size, Solution trans) throws JMException {
        double[][] mat1 = pop.getIndividual(task1, size);
        double[][] mat2 = pop.getIndividual(task2, size);
        DenoisingAutoencoder da = new DenoisingAutoencoder(mat1, mat2);
        da.train();
//        double[][] goodPops = Arrays.copyOfRange(mat1, 0, 25);
        double[][] res = da.predict(trans.getIndiVirables());
        return res;
    }

}
